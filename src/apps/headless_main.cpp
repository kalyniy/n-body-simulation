// headless_main.cpp
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>

#include "Simulation.h"
#include "PerformanceLogger.hpp"
#include "NaiveSimulation.h"
#include "BarnesHutSimulation.h"
#include "World.h"
#include "CheckpointManager.h"

#ifdef USE_MPI
#include <mpi.h>
#include "MpiNaiveSimulation.h"
#endif

#ifdef USE_CUDA
#include "CudaNaiveSimulation.h"
#ifdef USE_MPI
#include "CudaMpiNaiveSimulation.h"
#endif
#endif

static void printUsage(const char* prog) {
    std::cout
        << "Usage: " << prog << " [options]\n\n"
        << "Algorithm:\n"
        << "  --algo <naive|bh>       Algorithm (default: naive)\n"
        << "  --theta <float>         BH opening angle (default: 0.5)\n\n"
        << "Initial conditions:\n"
        << "  --random <N>    --solar    --hacc <dir>\n"
        << "  --plummer <N>   --galaxy <N>   --clusters <N>\n\n"
        << "Physics:\n"
        << "  --dt <f>  --G <f>  --softening <f>\n\n"
        << "Simulation:\n"
        << "  --steps <N>  --out <file>  --disable-output  --log-every <N>\n"
        << "  --help\n\n"
        << "Backend: "
#if defined(USE_CUDA) && defined(USE_MPI)
        << "CUDA+MPI"
#elif defined(USE_CUDA)
        << "CUDA"
#elif defined(USE_MPI)
        << "MPI"
#else
        << "Sequential"
#endif
        << "\n";
}

struct AlgoResult { std::unique_ptr<SimulationAlgorithm> algo; AlgorithmKind kind; };

static AlgoResult createAlgo(const std::string& name, float theta) {
    bool bh = (name=="bh"||name=="barnes-hut"||name=="barneshut");
    if (bh) {
        throw std::runtime_error("CUDA barns hut is not implemented yet");

    } else {
#if defined(USE_CUDA)&&defined(USE_MPI)
        std::cout<<"[algo] CudaNaive+MPI\n";
        return {std::make_unique<CudaMpiNaiveSimulation>(),AlgorithmKind::NaiveCudaMpi};
#elif defined(USE_CUDA)
        std::cout<<"[algo] CudaNaive\n";
        return {std::make_unique<CudaNaiveSimulation>(),AlgorithmKind::NaiveCuda};
#elif defined(USE_MPI)
        std::cout<<"[algo] NaiveMPI\n";
        return {std::make_unique<MpiNaiveSimulation>(),AlgorithmKind::NaiveMpi};
#else
        std::cout<<"[algo] Naive (seq)\n";
        return {std::make_unique<NaiveSimulation>(),AlgorithmKind::NaiveSeq};
#endif
    }
}

int main(int argc, char **argv)
{
#ifdef USE_MPI
    int rank, nproc;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
#else
    int rank=0, nproc=1;
#endif

    std::string algo_name="naive"; float theta=0.5f;
    std::string hacc_dir; size_t randomN=0, plummerN=0, galaxyN=0, clustersN=0;
    bool use_solar=false;
    float dt=0.5f, G=1.f, softening=2.f;
    size_t steps=1000, log_every=100; std::string out="output.txt";
    bool disable_output=false;

    for (int i=1;i<argc;++i) {
        std::string a=argv[i];
        if((a=="--algo"||a=="-a")&&i+1<argc) algo_name=argv[++i];
        else if(a=="--theta"&&i+1<argc) theta=std::stof(argv[++i]);
        else if(a=="--hacc"&&i+1<argc){hacc_dir=argv[++i];use_solar=false;}
        else if(a=="--random"&&i+1<argc){randomN=std::stoul(argv[++i]);use_solar=false;}
        else if(a=="--plummer"&&i+1<argc){plummerN=std::stoul(argv[++i]);use_solar=false;}
        else if(a=="--galaxy"&&i+1<argc){galaxyN=std::stoul(argv[++i]);use_solar=false;}
        else if(a=="--clusters"&&i+1<argc){clustersN=std::stoul(argv[++i]);use_solar=false;}
        else if(a=="--solar") use_solar=true;
        else if(a=="--dt"&&i+1<argc) dt=std::stof(argv[++i]);
        else if(a=="--G"&&i+1<argc) G=std::stof(argv[++i]);
        else if(a=="--softening"&&i+1<argc) softening=std::stof(argv[++i]);
        else if(a=="--steps"&&i+1<argc) steps=std::stoul(argv[++i]);
        else if(a=="--out"&&i+1<argc) out=argv[++i];
        else if(a=="--disable-output") disable_output=true;
        else if(a=="--log-every"&&i+1<argc) log_every=std::stoul(argv[++i]);
        else if(a=="--help"||a=="-h"){if(rank==0)printUsage(argv[0]);
#ifdef USE_MPI
            MPI_Finalize();
#endif
            return 0;}
    }

    SimParams params; params.G=G; params.dt=dt; params.min_r2=softening;
    auto [algorithm, algo_kind] = createAlgo(algo_name, theta);
    NBodySimulation sim(std::move(algorithm), params);

    size_t n_global = 0;

#ifdef USE_MPI
    bool domain_decomposed =
        (algo_kind==AlgorithmKind::NaiveMpi||algo_kind==AlgorithmKind::BarnesHutMpi||
         algo_kind==AlgorithmKind::NaiveCudaMpi||algo_kind==AlgorithmKind::BarnesHutCudaMpi);
#endif

    // Generate particles
    auto gen = [&](auto fn) {
#ifdef USE_MPI
        if(rank==0){fn();n_global=sim.particles().size();std::printf("Rank 0: %zu particles\n",n_global);}
#else
        fn(); n_global=sim.particles().size();
#endif
    };

    if(!hacc_dir.empty()) gen([&]{sim.loadHACC(hacc_dir);});
    else if(plummerN>0) gen([&]{sim.generatePlummerSphere(plummerN,100.f,float(plummerN));});
    else if(galaxyN>0) gen([&]{sim.generateGalaxyDisk(galaxyN,300.f,50.f);});
    else if(clustersN>0) gen([&]{sim.generateClusters(3,clustersN,500.f);});
    else if(randomN>0) gen([&]{sim.generateRandom(randomN,WORLD_WIDTH,WORLD_HEIGHT,WORLD_DEPTH);});
    else if(use_solar) gen([&]{sim.setupSolarSystem(WORLD_WIDTH,WORLD_HEIGHT,WORLD_DEPTH);});
    else {
        if(rank==0) std::cerr<<"No initial condition. Use --help\n";
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
        return 1;
    }

    // MPI distribution
#ifdef USE_MPI
    MPI_Bcast(&n_global,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
    // For MPI we need to scatter SoA arrays
    std::vector<int> counts(nproc), displs(nproc,0);
    if(domain_decomposed){
        int base=n_global/nproc, rem=n_global%nproc;
        for(int r=0;r<nproc;++r) counts[r]=base+(r<rem?1:0);
        for(int r=1;r<nproc;++r) displs[r]=displs[r-1]+counts[r-1];
        int n_local=counts[rank];

        // Scatter each SoA array separately
        ParticleSoA global_copy;
        if(rank==0) global_copy = sim.particles();

        sim.particles().resize(n_local);
        auto scatter_arr = [&](std::vector<float>& local, const std::vector<float>& global) {
            MPI_Scatterv(rank==0?global.data():nullptr,counts.data(),displs.data(),MPI_FLOAT,
                         local.data(),n_local,MPI_FLOAT,0,MPI_COMM_WORLD);
        };
        scatter_arr(sim.particles().pos_x, global_copy.pos_x);
        scatter_arr(sim.particles().pos_y, global_copy.pos_y);
        scatter_arr(sim.particles().pos_z, global_copy.pos_z);
        scatter_arr(sim.particles().vel_x, global_copy.vel_x);
        scatter_arr(sim.particles().vel_y, global_copy.vel_y);
        scatter_arr(sim.particles().vel_z, global_copy.vel_z);
        scatter_arr(sim.particles().acc_x, global_copy.acc_x);
        scatter_arr(sim.particles().acc_y, global_copy.acc_y);
        scatter_arr(sim.particles().acc_z, global_copy.acc_z);
        scatter_arr(sim.particles().mass, global_copy.mass);
        std::printf("Rank %d: local n = %d\n",rank,n_local);
    } else {
        // Broadcast full arrays
        if(rank!=0) sim.particles().resize(n_global);
        auto bcast_arr = [&](std::vector<float>& v) {
            MPI_Bcast(v.data(),n_global,MPI_FLOAT,0,MPI_COMM_WORLD);
        };
        bcast_arr(sim.particles().pos_x); bcast_arr(sim.particles().pos_y); bcast_arr(sim.particles().pos_z);
        bcast_arr(sim.particles().vel_x); bcast_arr(sim.particles().vel_y); bcast_arr(sim.particles().vel_z);
        bcast_arr(sim.particles().mass);
        sim.particles().zeroAccelerations();
    }
#else
    if(n_global==0){std::cerr<<"No particles\n";return 1;}
#endif

    // Checkpoint setup
    PerformanceLogger logger(out, RunMode::Run);
    if(rank==0&&!logger.ok()) std::cerr<<"Warning: can't open "<<out<<"\n";

    if(!disable_output){
#ifdef USE_MPI
        if(rank==0)
#endif
        {
            auto* cp=CheckpointManager::getInstance();
            cp->setFilePath("simulation_output.bin");
            SimulationOutputHeader h; h.n_particles=n_global; h.target_steps=steps; h.passed_steps=0;
            cp->write_header(h);
#ifdef USE_MPI
            // Need global particles for masses — but we scattered. Re-gather or use saved copy.
            // For simplicity, rank 0 still has the global copy from generation
            ParticleSoA global_copy; global_copy.resize(n_global);
            // Gather masses
            MPI_Gatherv(sim.particles().mass.data(),counts[rank],MPI_FLOAT,
                        global_copy.mass.data(),counts.data(),displs.data(),MPI_FLOAT,0,MPI_COMM_WORLD);
            cp->write_masses(global_copy);
#else
            cp->write_masses(sim.particles());
#endif
        }
#ifdef USE_MPI
        else {
            // Non-zero ranks still need to participate in Gatherv
            MPI_Gatherv(sim.particles().mass.data(),counts[rank],MPI_FLOAT,
                        nullptr,nullptr,nullptr,MPI_FLOAT,0,MPI_COMM_WORLD);
        }
#endif
    }

    if(rank==0) {
        std::cout<<"\n=== Run Configuration ===\n"
                 <<"  Algorithm:   "<<algo_name<<"\n"
                 <<"  Particles:   "<<n_global<<"\n"
                 <<"  Steps:       "<<steps<<"\n"
                 <<"  dt:          "<<params.dt<<"\n"
                 <<"  G:           "<<params.G<<"\n"
                 <<"  Softening^2: "<<params.min_r2<<"\n"
                 <<"  Processes:   "<<nproc<<"\n"
#ifdef USE_CUDA
                 <<"  CUDA:        enabled\n"
#else
                 <<"  CUDA:        disabled\n"
#endif
                 <<"=========================\n\n";
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    auto t0 = std::chrono::high_resolution_clock::now();

#ifdef USE_MPI
    int n_local = sim.particles().size();
    ParticleSoA global_buf;
    if(rank==0&&domain_decomposed) global_buf.resize(n_global);
#endif

    for(size_t step=0;step<steps;++step) {
        sim.step();

#ifdef USE_MPI
        if(domain_decomposed) {
            // Gather SoA arrays for checkpoint
            if(!disable_output || step%log_every==0) {
                auto gather=[&](const std::vector<float>& local, std::vector<float>& global) {
                    MPI_Gatherv(local.data(),n_local,MPI_FLOAT,
                                rank==0?global.data():nullptr,counts.data(),displs.data(),MPI_FLOAT,0,MPI_COMM_WORLD);
                };
                gather(sim.particles().pos_x,global_buf.pos_x);
                gather(sim.particles().pos_y,global_buf.pos_y);
                gather(sim.particles().pos_z,global_buf.pos_z);
            }
            if(rank==0) {
                if(!disable_output) CheckpointManager::getInstance()->write_step(global_buf);
                if(step%log_every==0) std::cout<<"Step "<<step<<" / "<<steps<<"\n";
            }
        } else {
            if(rank==0) {
                if(!disable_output) CheckpointManager::getInstance()->write_step(sim.particles());
                if(step%log_every==0) std::cout<<"Step "<<step<<" / "<<steps<<"\n";
            }
        }
#else
        if(!disable_output) CheckpointManager::getInstance()->write_step(sim.particles());
        if(step%log_every==0) std::cout<<"Step "<<step<<" / "<<steps<<"\n";
#endif
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1-t0).count();
    if(rank==0) {
        logger.log(n_global, steps, elapsed, nproc);
        std::cout<<"Simulation completed in "<<elapsed<<" seconds\n";
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
