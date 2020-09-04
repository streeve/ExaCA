#include "header.h"
using namespace std;

void RunProgram_Reduced(int id, int np, int ierr, string InputFile) {
    
    
    double InitTime = MPI_Wtime();
    
    int nx, ny, nz, DecompositionStrategy, NumberOfLayers, LayerHeight, TempFilesInSeries, NumberOfTruchasRanks;
    bool FilesToPrint[4];
    double HT_deltax, deltax, deltat, FractSurfaceSitesActive, G, R, AConst, BConst, CConst, DConst, FreezingRange, NMax, dTN, dTsigma;
    string SubstrateFileName, tempfile, SimulationType, OutputFile, GrainOrientationFile, TemperatureDataSource, BurstBuffer, ExtraWalls, PathToOutput;
    
    // Read input data
    InputReadFromFile(id, InputFile, SimulationType, DecompositionStrategy, AConst, BConst, CConst, DConst, FreezingRange, deltax, NMax, dTN, dTsigma, OutputFile, GrainOrientationFile, tempfile, TempFilesInSeries, BurstBuffer, ExtraWalls, HT_deltax, TemperatureDataSource, deltat, NumberOfLayers, LayerHeight, SubstrateFileName, G, R, nx, ny, nz, FractSurfaceSitesActive,PathToOutput,NumberOfTruchasRanks,FilesToPrint);
    
    // Grid decomposition
    int ProcessorsInXDirection, ProcessorsInYDirection;
    // Variables characterizing local processor grids relative to global domain
    int MyXSlices, MyXOffset, MyYSlices, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut;
    // Neighbor lists for cells
    int NeighborX[26], NeighborY[26], NeighborZ[26], ItList[9][26];
    float XMin, YMin, ZMin, XMax, YMax, ZMax; // OpenFOAM simulation bounds (if using OpenFOAM data)
    float* ZMinLayer = new float[NumberOfLayers];
    float* ZMaxLayer = new float[NumberOfLayers];
    int* FinishTimeStep = new int[NumberOfLayers];
    
    // Temporary data structure for storing temperature data from file(s)
    // Initial estimate for size
    vector <float> RawData(1000000);
    // Contains "NumberOfLayers" values corresponding to the location within "RawData" of the first data element in each layer
    int* FirstValue = new int[NumberOfLayers];

    // Initialization of the grid and decomposition, along with deltax and deltat
    // Read in temperature data
    ParallelMeshInit(DecompositionStrategy, NeighborX, NeighborY, NeighborZ, ItList, SimulationType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, TemperatureDataSource, LayerHeight, NumberOfLayers, TempFilesInSeries, ZMinLayer, ZMaxLayer, FirstValue, RawData, BurstBuffer,NumberOfTruchasRanks);
    
    long int LocalDomainSize = MyXSlices*MyYSlices*nz; // Number of cells on this MPI rank

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Mesh initialized" << endl;
    
    // Temperature fields characterized by these variables:
     ViewI CritTimeStep_G("CritTimeStep_G", LocalDomainSize);
     ViewI LayerID_G("LayerID_G", LocalDomainSize);
     ViewF UndercoolingChange_G("UndercoolingChange_G",LocalDomainSize);
     ViewF UndercoolingCurrent_G("UndercoolingCurrent_G",LocalDomainSize);
     ViewI::HostMirror CritTimeStep_H = Kokkos::create_mirror_view( CritTimeStep_G );
     ViewI::HostMirror LayerID_H = Kokkos::create_mirror_view( LayerID_G );
     ViewF::HostMirror UndercoolingChange_H = Kokkos::create_mirror_view( UndercoolingChange_G );
     ViewF::HostMirror UndercoolingCurrent_H = Kokkos::create_mirror_view( UndercoolingCurrent_G );
     bool* Melted = new bool[LocalDomainSize];

    // By default, the active domain bounds are the same as the global domain bounds
    // For multilayer problems, this is not the case and ZBound_High and nzActive will be adjusted in TempInit to account for only the first layer of solidification
    int ZBound_Low;
    int ZBound_High;
    int nzActive;
    
    // Initialize the temperature fields
    TempInit(-1, TempFilesInSeries, G, R, DecompositionStrategy,NeighborX, NeighborY, NeighborZ, ItList, SimulationType, ierr, id, np, MyXSlices, MyYSlices, MyXOffset, MyYOffset, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyLeftOut, MyRightIn, MyRightOut, deltax, HT_deltax, deltat, nx, ny, nz, ProcessorsInXDirection, ProcessorsInYDirection,  CritTimeStep_H, UndercoolingChange_H, UndercoolingCurrent_H, tempfile, XMin, XMax, YMin, YMax, ZMin, ZMax, Melted, TemperatureDataSource, ZMinLayer, ZMaxLayer, LayerHeight, NumberOfLayers, nzActive, ZBound_Low, ZBound_High, FinishTimeStep, FreezingRange, LayerID_H, FirstValue, RawData, BurstBuffer);
    // Delete temporary data structure for temperature data read
    RawData.clear();
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with temperature field initialization, active domain size is " << nzActive << " out of " << nz << " cells in the Z direction" << endl;

    int LocalActiveDomainSize = MyXSlices*MyYSlices*nzActive; // Number of active cells on this MPI rank
    
    // PrintTempValues(id,np,nx,ny,nz, MyXSlices, MyYSlices, ProcessorsInXDirection, ProcessorsInYDirection, CritTimeStep_H, UndercoolingChange_H, DecompositionStrategy,PathToOutput);
    
    
    int NGrainOrientations = 10000; // Number of grain orientations considered in the simulation
    
    
    float* GrainUnitVector = new float[9*NGrainOrientations];
    int* GrainOrientation = new int[NGrainOrientations];
    
    // Initialize grain orientations
    OrientationInit(id, NGrainOrientations, GrainOrientation, GrainUnitVector, GrainOrientationFile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with orientation initialization " << endl;
    
    // CA cell variables
    ViewI GrainID_G("GrainID_G",LocalDomainSize);
    ViewI CellType_G("CellType_G",LocalDomainSize);
    ViewI::HostMirror GrainID_H = Kokkos::create_mirror_view( GrainID_G );
    ViewI::HostMirror CellType_H = Kokkos::create_mirror_view( CellType_G );

    // Variables characterizing the active cell region within each rank's grid
    ViewF DiagonalLength_G("DiagonalLength_G",LocalActiveDomainSize);
    ViewF CritDiagonalLength_G("CritDiagonalLength_G",26*LocalActiveDomainSize);
    ViewF DOCenter_G("DOCenter_G",3*LocalActiveDomainSize);
    ViewF::HostMirror DiagonalLength_H = Kokkos::create_mirror_view( DiagonalLength_G );
    ViewF::HostMirror CritDiagonalLength_H = Kokkos::create_mirror_view( CritDiagonalLength_G );
    ViewF::HostMirror DOCenter_H = Kokkos::create_mirror_view( DOCenter_G );

    // Initialize the grain structure
    int PossibleNuclei_ThisRank, NextLayer_FirstNucleatedGrainID;
    
    GrainInit(-1, LayerHeight, SimulationType, SubstrateFileName, FractSurfaceSitesActive, NGrainOrientations, DecompositionStrategy, ProcessorsInXDirection, ProcessorsInYDirection, nx, ny, nz, MyXSlices, MyYSlices, MyXOffset, MyYOffset, id, np, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, ItList, NeighborX, NeighborY, NeighborZ, GrainOrientation, GrainUnitVector, DiagonalLength_H, CellType_H, GrainID_H, CritDiagonalLength_H, DOCenter_H, CritTimeStep_H, UndercoolingChange_H, Melted, deltax, NMax, NextLayer_FirstNucleatedGrainID, PossibleNuclei_ThisRank, ZBound_High, NumberOfLayers, TempFilesInSeries, HT_deltax, deltat, XMin, XMax, YMin, YMax, ZMin, ZMax,tempfile,TemperatureDataSource, ZBound_Low, ExtraWalls);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Grain struct initialized" << endl;

    ViewI NucleationTimes_G("NucleationTimes_G",PossibleNuclei_ThisRank);
    ViewI NucleiLocation_G("NucleiLocation_G",PossibleNuclei_ThisRank);
    ViewI::HostMirror NucleationTimes_H = Kokkos::create_mirror_view( NucleationTimes_G );
    ViewI::HostMirror NucleiLocation_H = Kokkos::create_mirror_view( NucleiLocation_G );
    
    // Update nuclei on ghost nodes, fill in nucleation data structures, and assign nucleation undercooling values to potential nucleation events
    NucleiInit(DecompositionStrategy, MyXSlices, MyYSlices, nz, id, dTN, dTsigma, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, PossibleNuclei_ThisRank, NucleiLocation_H, NucleationTimes_H, GrainOrientation, CellType_H, GrainID_H, CritTimeStep_H, UndercoolingChange_H);

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Done with nuclei initialization " << endl;
    
    // Normalize solidification parameters
    AConst = AConst*deltat/deltax;
    BConst = BConst*deltat/deltax;
    CConst = CConst*deltat/deltax;
    int cycle;
    
    // Buffers for ghost node data (fixed size)
    int BufSizeX, BufSizeY, BufSizeZ;
    if (DecompositionStrategy == 1) {
        BufSizeX = MyXSlices;
        BufSizeY = 0;
        BufSizeZ = nzActive;
    }
    else {
        BufSizeX = MyXSlices-2;
        BufSizeY = MyYSlices-2;
        BufSizeZ = nzActive;
    }
    Buffer2D BufferA("BufferA",BufSizeX*BufSizeZ,5);
    Buffer2D BufferB("BufferB",BufSizeX*BufSizeZ,5);
    Buffer2D BufferC("BufferC",BufSizeY*BufSizeZ,5);
    Buffer2D BufferD("BufferD",BufSizeY*BufSizeZ,5);
    Buffer2D BufferE("BufferE",BufSizeZ,5);
    Buffer2D BufferF("BufferF",BufSizeZ,5);
    Buffer2D BufferG("BufferG",BufSizeZ,5);
    Buffer2D BufferH("BufferH",BufSizeZ,5);
    Buffer2D BufferAR("BufferAR",BufSizeX*BufSizeZ,5);
    Buffer2D BufferBR("BufferBR",BufSizeX*BufSizeZ,5);
    Buffer2D BufferCR("BufferCR",BufSizeY*BufSizeZ,5);
    Buffer2D BufferDR("BufferDR",BufSizeY*BufSizeZ,5);
    Buffer2D BufferER("BufferER",BufSizeZ,5);
    Buffer2D BufferFR("BufferFR",BufSizeZ,5);
    Buffer2D BufferGR("BufferGR",BufSizeZ,5);
    Buffer2D BufferHR("BufferHR",BufSizeZ,5);
    
    // Copy view data to GPU
    Kokkos::deep_copy( GrainID_G, GrainID_H );
    Kokkos::deep_copy( CellType_G, CellType_H );
    Kokkos::deep_copy( DiagonalLength_G, DiagonalLength_H );
    Kokkos::deep_copy( CritDiagonalLength_G, CritDiagonalLength_H );
    Kokkos::deep_copy( DOCenter_G, DOCenter_H );
    Kokkos::deep_copy( CritTimeStep_G, CritTimeStep_H );
    Kokkos::deep_copy( LayerID_G, LayerID_H );
    Kokkos::deep_copy( UndercoolingChange_G, UndercoolingChange_H );
    Kokkos::deep_copy( UndercoolingCurrent_G, UndercoolingCurrent_H );
    Kokkos::deep_copy( NucleiLocation_G, NucleiLocation_H );
    Kokkos::deep_copy( NucleationTimes_G, NucleationTimes_H );

    
    // Locks for cell capture
    // 0 = cannot be captured, 1 = can be capured
    ViewI Locks("Locks",LocalActiveDomainSize);
    Kokkos::parallel_for("LockInit",LocalActiveDomainSize, KOKKOS_LAMBDA (const int& D3D1ConvPosition) {
        int RankZ = floor(D3D1ConvPosition/(MyXSlices*MyYSlices));
        int Rem = D3D1ConvPosition % (MyXSlices*MyYSlices);
        int RankX = floor(Rem/MyYSlices);
        int RankY = Rem % MyYSlices;
        int GlobalZ = ZBound_Low + RankZ;
        int GlobalD3D1ConvPosition = GlobalZ*MyXSlices*MyYSlices + RankX*MyYSlices + RankY;
        if ((CellType_G(GlobalD3D1ConvPosition) == Delayed)||(CellType_G(GlobalD3D1ConvPosition) == LiqSol)||(CellType_G(GlobalD3D1ConvPosition) == Liquid)) Locks(D3D1ConvPosition) = 1;
        else Locks(D3D1ConvPosition) = 0;
    });
    
    if (np > 1) {
        // Ghost nodes for initial microstructure state
        GhostNodesInit_GPU(DecompositionStrategy, MyXSlices, MyYSlices, GrainID_G, CellType_G, DOCenter_G, DiagonalLength_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, LocalActiveDomainSize, ZBound_Low);
        if (DecompositionStrategy == 1) GhostNodes1D_GPU(0, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
        else GhostNodes2D_GPU(0, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
    }

    double InitTime2 = MPI_Wtime();
    if (id == 0) cout << "Data initialized: Time spent: " << InitTime2-InitTime << " s" << endl;

    cycle = 0;
    
    for (int layernumber=0; layernumber<NumberOfLayers; layernumber++) {
        
        int nn = 0; // Counter for the number of nucleation events
        int XSwitch = 0;
        double LayerTime1 = MPI_Wtime();
        
        // Loop continues until all liquid cells claimed by solid grains
        do {
            cycle++;

            
            //if ((layernumber == 0)&&(id == 0)) cout << " CYCLE " << cycle << endl;
            // Update cells on GPU - undercooling and diagonal length updates, nucleation
           Nucleation(id, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, cycle, nn, CritTimeStep_G, CellType_G, UndercoolingCurrent_G, UndercoolingChange_G, NucleiLocation_G, NucleationTimes_G, GrainID_G, GrainOrientation, DOCenter_G, NeighborX,  NeighborY, NeighborZ, GrainUnitVector, CritDiagonalLength_G, DiagonalLength_G, NGrainOrientations, PossibleNuclei_ThisRank, Locks, ZBound_Low, layernumber, LayerID_G);

            // Update cells on GPU - new active cells, solidification of old active cells
            CellCapture(id, np, cycle, DecompositionStrategy, LocalActiveDomainSize, MyXSlices, MyYSlices, nz, AConst, BConst, CConst, DConst, MyXOffset, MyYOffset, ItList, NeighborX, NeighborY, NeighborZ, CritTimeStep_G, UndercoolingCurrent_G, UndercoolingChange_G,  GrainUnitVector, CritDiagonalLength_G, DiagonalLength_G, GrainOrientation, CellType_G, DOCenter_G, GrainID_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, Locks, ZBound_Low, nzActive, layernumber, LayerID_G);
            //cout << "ID = " << id << " waiting" << endl;
           //  MPI_Barrier(MPI_COMM_WORLD);
           // if ((layernumber == 0)&&(id == 0)) cout << " CYCLE " << cycle << endl;
            if (np > 1) {
            // Update ghost nodes
                if (DecompositionStrategy == 1) GhostNodes1D_GPU(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
                else GhostNodes2D_GPU(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
            }
            //if ((layernumber == 1)&&(id == 0)) cout << " CYCLE " << cycle << endl;
            if (cycle % 1000 == 0) {
                IntermediateOutputAndCheck(id, cycle, MyXSlices, MyYSlices, LocalDomainSize, LocalActiveDomainSize, nn, XSwitch, CellType_G, CritTimeStep_G, SimulationType, FinishTimeStep, layernumber, NumberOfLayers, ZBound_Low, Locks, LayerID_G);
            }
            

            //if (cycle == 20000) XSwitch = 1;
            //if (cycle == 20000)     PrintCT( id,  np,  nx,  ny,  nz,  MyXSlices,  MyYSlices,  ProcessorsInXDirection,  ProcessorsInYDirection,  CellType_H,  OutputFile,  DecompositionStrategy);
            
        } while(XSwitch == 0);


        if (layernumber != NumberOfLayers-1) {
             // Determine new active cell domain size and offset from bottom of global domain
             int ZShift;
             DomainShiftAndResize(id, MyXSlices, MyYSlices, ZShift, ZBound_Low, ZBound_High, nzActive, LocalDomainSize, LocalActiveDomainSize, BufSizeZ, LayerHeight, CellType_G, layernumber, LayerID_G);

             // Resize active cell data structures
             Kokkos::resize(DiagonalLength_G,LocalActiveDomainSize);
             Kokkos::resize(DOCenter_G,3*LocalActiveDomainSize);
             Kokkos::resize(CritDiagonalLength_G,26*LocalActiveDomainSize);
             Kokkos::resize(Locks,LocalActiveDomainSize);
            
             Kokkos::resize(BufferA,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferB,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferC,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferD,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferE,BufSizeZ,5);
             Kokkos::resize(BufferF,BufSizeZ,5);
             Kokkos::resize(BufferG,BufSizeZ,5);
             Kokkos::resize(BufferH,BufSizeZ,5);
            
             Kokkos::resize(BufferAR,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferBR,BufSizeX*BufSizeZ,5);
             Kokkos::resize(BufferCR,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferDR,BufSizeY*BufSizeZ,5);
             Kokkos::resize(BufferER,BufSizeZ,5);
             Kokkos::resize(BufferFR,BufSizeZ,5);
             Kokkos::resize(BufferGR,BufSizeZ,5);
             Kokkos::resize(BufferHR,BufSizeZ,5);
            
             MPI_Barrier(MPI_COMM_WORLD);
             if (id == 0) cout << "Resize executed" << endl;
            
             // Update active cell data structures for simulation of next layer
             LayerSetup(SubstrateFileName, layernumber, LayerHeight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, LocalDomainSize, LocalActiveDomainSize, GrainOrientation, NGrainOrientations, GrainUnitVector, NeighborX, NeighborY, NeighborZ, id, np, DiagonalLength_G, CellType_G, GrainID_G, CritDiagonalLength_G, DOCenter_G, CritTimeStep_G, UndercoolingChange_G, UndercoolingCurrent_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, TempFilesInSeries, ZBound_Low, ZBound_High, ZMin, deltax, nzActive, Locks);

            if (id == 0) cout << "New layer setup, GN dimensions are " << BufSizeX << " " << BufSizeY << " " << BufSizeZ << endl;
            // Update ghost nodes for grain locations and attributes
            GhostNodesInit_GPU(DecompositionStrategy, MyXSlices, MyYSlices, GrainID_G, CellType_G, DOCenter_G, DiagonalLength_G, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufSizeX, BufSizeY, LocalActiveDomainSize, ZBound_Low);
            MPI_Barrier(MPI_COMM_WORLD);
            if (id == 0) cout << "New layer ghost nodes initialized" << endl;
             if (np > 1) {
             // Update ghost nodes
                 if (DecompositionStrategy == 1) GhostNodes1D_GPU(cycle, id, MyLeft, MyRight, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G,GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferAR, BufferBR, BufSizeX,  BufSizeY, BufSizeZ, Locks, ZBound_Low);
                 else GhostNodes2D_GPU(cycle, id, MyLeft, MyRight, MyIn, MyOut, MyLeftIn, MyRightIn, MyLeftOut, MyRightOut, MyXSlices, MyYSlices, MyXOffset, MyYOffset, nz, NeighborX, NeighborY, NeighborZ, CellType_G, DOCenter_G, GrainID_G, GrainUnitVector, GrainOrientation, DiagonalLength_G, CritDiagonalLength_G, NGrainOrientations, BufferA, BufferB, BufferC, BufferD, BufferE, BufferF, BufferG, BufferH, BufferAR, BufferBR, BufferCR, BufferDR, BufferER, BufferFR, BufferGR, BufferHR, BufSizeX, BufSizeY, BufSizeZ, Locks, ZBound_Low);
             }
             XSwitch = 0;
             MPI_Barrier(MPI_COMM_WORLD);
             double LayerTime2 = MPI_Wtime();
             cycle = 0;
             if (id == 0) cout << "Time for layer number " << layernumber << " was " << LayerTime2-LayerTime1 << " s, starting layer " << layernumber+1 << endl;
         }
         else {
            MPI_Barrier(MPI_COMM_WORLD);
            double LayerTime2 = MPI_Wtime();
            if (id == 0) cout << "Time for final layer was " << LayerTime2-LayerTime1 << " s" << endl;
         }

    }

    double InitTime3 = MPI_Wtime();

    // Copy GPU results for GrainID back to CPU for printing to file(s)
    Kokkos::deep_copy( GrainID_H, GrainID_G );
    Kokkos::deep_copy( CellType_H, CellType_G );

    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) cout << "Collecting data on rank 0 and printing to files" << endl;
    CollectGrainData(id,np,nx,ny,nz,MyXSlices,MyYSlices, MyXOffset, MyYOffset, ProcessorsInXDirection, ProcessorsInYDirection, GrainID_H,GrainOrientation,GrainUnitVector,OutputFile,DecompositionStrategy,NGrainOrientations,Melted,PathToOutput,FilesToPrint,deltax);
    
    double InitTime4 = MPI_Wtime();
    if (id == 0) {
        cout << "===================================================================================" << endl;
        cout << "Having run with = " << np << " processors" << endl;
        cout << "Output written at cycle = " << cycle << endl;
        cout << "Total time = " << InitTime4 - InitTime << endl;
        cout << "Time spent initializing data = " << InitTime2-InitTime << " s" << endl;
        cout << "Time spent performing CA calculations = " << InitTime3-InitTime2 << " s" << endl;
        cout << "Time spent collecting and printing output data = " << InitTime4-InitTime3 << " s" << endl;
        cout << "===================================================================================" << endl;
    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "ID = " << id << " ready to exit" << endl;
//    double GlobalT;
//    MPI_Reduce(&CoolTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in cooling = " << GlobalT << " s" << endl;
//    MPI_Reduce(&NucTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in nucleation = " << GlobalT << " s" << endl;
//    MPI_Reduce(&GrowTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in growth = " << GlobalT << " s" << endl;
//    MPI_Reduce(&CaptTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in capture = " << GlobalT << " s" << endl;
//    MPI_Reduce(&DeactTime,&GlobalT,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    if (id == 0) cout << "Time spent in deactivation = " << GlobalT << " s" << endl;
    
}
