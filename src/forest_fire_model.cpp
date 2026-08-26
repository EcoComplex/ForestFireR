/*********************************************************************/
//  Spatial forest-fire model
/*********************************************************************/
// Stochastic model based on a cellular automaton of forest-fire dynamics.
// The model is implemented by the Gillespie algorithm on a 2D grid
// The model has 5 possible states:
// **S: Possible States, 
//      S=0 -> Empty
//      S=1 -> Vegetation 1
//      S=2 -> Vegetation 2
//      S=3 -> Fire
//      S=4 -> Empty_Vegetation
//      S=5 -> Empty_Fire
// The model computes the evolution of the system from an initial condition 
// and a set of transition probabilities.
//*********************************************************************/
// Model parameters:
// L: Grid size, N=LxL
// PeriodicBoundaryConditions={true,false}: Periodic boundary conditions
//
//T: time (years)
//q22: parameter p indicating the heterogeneity in the initial distribution of species
//density2: initial density of the invader species
//
// Transition rates
// L_01
// L_02
// L_20
// L_21
// L_12
// L_30
// Lsp_13
// Lsp_23
// Lig_13
// Lig_23
// Lrg_02
// xim
// etam
//
// Output parameters:
// time_sim: time in the simulation
// n0,n1,n2,n3: Species densities
//
//*********************************************************************/
//---------------------Functions-----------------------------------
//
// -> InitialConditionNonHomogeneous(density2,q22); 
// Assign an initial condition for the spatial 
// configuration S[][] according to the initial 
// density of the invader species and the parameter p
//------Code Example (1):
// q22=1; 
// density2=0.1;
// InitialConditionNonHomogeneous(density2,q22);
// for(int i=0;i<L;i++)
// {
//     for(int j=0;j<L;j++)
//     {
//         cout<<i<<" "<<j<<" "<<S[i][j]<<endl;
//     }
// } 
//------End Code Example (1):  
//
// -> SpatialModelSimplified(T);  
// Run the model from a spatial initial condition S[][],
// fixed transition probabilities, and a time T
//------Code Example (2):
// q22=1; 
// density2=0.1;
// InitialConditionNonHomogeneous(density2,q22);
// T=200
// SpatialModelSimplified(T);
// cout<<time_sim<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl; 
//------End Code Example (2):        
// 
// -> MeanField(T);
// -> MeanField_Guillespie(T);
// Compute the mean-field and the the stochastic version
// of the mean-field by the gillespie algorithm
/*********************************************************************/
//  Computed by Rebeca de la Fuente :  7/11/2024      
//  science.rdelafuente@hotmail.com                           
/*********************************************************************/



#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <Rcpp.h>

using namespace std;

// #include "Struct.h"
// #include "Function.h"
// #include "parallel.h"

struct Vector {
    long double x, y, z, q, s, l, t;
};



//System Size
int L=200; 
int N=L*L; 
long double system_size=N;

//Model States
int **S = nullptr;

//Model Parameters
long double L_10=0;
long double L_01=0.03;
long double L_20=0;
long double L_02=0.03;
long double L_21=0.01;
long double L_12=0.005; 
long double L_30=1000000;
long double Lsp_13=1000000;
long double Lig_13=0; 
long double Lig_23=0.0001;

long double xim,etam;
long double Lsp_23=0;
long double Lrg_02=0;

long double Lrg_01=0;
long double Lr_12=0; 
long double Lr_21=0; 
long double Lr_01=0;
long double Lr_02=0;

//Particle Densities
long double N0,N1,N2,N3,N4,N5;
long double n0,n1,n2,n3,n4,n5;



//Reactions
int r=13; //total number of reactions
//-----Vegetation Reactions
//0 -> Empty_v to V1
//1 -> V1 to Empty_v
//2 -> Empty_v to V2
//3 -> V2 to Empty_v
//4 -> V1 to V2lag time
//5 -> V2 to V1
//--------Fire Reactions
//6 -> V1 to F
//7 -> V2 to F
//8 -> F to Emtpy_f
//9 -> Empty_f to V1
//10 -> Empty_f to V2
//11 -> Empty_f to V1 (vegetation dynamics)
//12 -> Empty_f to V2 (vegetation dynamics)


//Transition rates
long double **A = nullptr;


//Functions
//ReverseGridMap(ix,jx,L,k);
//GridMap(ix,jx,L,k);
void VegetationFireModel_Rates(int &);
void MeanField_Guillespie(long double &);
Vector ForestFire(Vector &);
Vector zRK4(Vector &);
void MeanField(long double &);
void SpatialModel(long double &);
void SpatialModelSimplified(long double &);
void SpatialModel_Movie(long double &);
long double neighbours(int &,int &,int &);
long double full_neighbours(int &,int &,int &);
vector<int> Cluster(int &);
long double neighbours_PeriodicBoundaryConditions(int &,int &,int &);
void InitialCondition(long double &);
void InitialConditionInClusters(long double &);
bool Boundaries(int &, int &);
void InitialConditionFromRandomness(long double &, int &);
long double Entropy(long double &);
// state1_param/state2_param generalize the original hardcoded "background
// = 1, pattern = 2" landscape (native vs. invader) to an arbitrary pair of
// grid states, and reset_grid controls whether the whole grid is first
// filled with state1_param (the original, single-layer behavior) or left
// as-is so a pattern can be grown into the surviving background of a prior
// call -- see FillGrid()/TallyGridToCompartments() below and the
// generate_landscape_layers_cpp/simulate_spatial_from_grid_cpp Rcpp
// exports, added to build multi-species / fire-over-invader landscapes
// (Figs. 8-11) that a single background/pattern pair can't express.
void InitialConditionNonHomogeneous(long double &, long double &, int state1_param = 1, int state2_param = 2, bool reset_grid = true);
long double H();

//---- Added for R/Rcpp integration ----
void GridMap(int ix, int jx, int Lp, int &k);
void ReverseGridMap(int &ix, int &jx, int Lp, int k);
void AllocateGrids(int newL);
void SeedRNG(int seed);
// check_extinction reproduces the original model's early stop once native
// density n1 drops below 0.0001 (a sensible default when n1 starts high,
// as in every figure the paper itself reproduces this way). It's wrong,
// though, for an initial condition built with no native vegetation at all
// -- e.g. Figs. 10-11's fire-over-invader landscape -- where n1 is 0 from
// the very first step and the run would otherwise stop after a single
// event. simulate_spatial_cpp always passes true (unchanged behavior);
// simulate_spatial_from_grid_cpp lets the caller turn it off.
void SpatialModelSimplified_Rec(long double &T, long double &record_dt, bool check_extinction = true);
void MeanField_Rec(long double &T, long double &dt, double n1_0, double n2_0, double n3_0, double n4_0, double n5_0);
void MeanFieldGillespie_Rec(long double &T, long double &record_dt, long double &sysN);
void FillGrid(int state);
void TallyGridToCompartments();
//---------------------------------------


long double xc;
long double parameter;
long double xvacio;
bool PeriodicBoundaryConditions;
long double time_sim;









int forest_fire_model_unused_main() 
{
return 0;
#if 0

S = new int *[L];
for (int i = 0; i < L; i++){S[i] = new int [L];}
A = new long double *[N];
for (int i = 0; i < N; i++){A[i] = new long double [r];}

//Providing a seed value
srand((unsigned) time(NULL));

//Periodic Boundary Conditions
PeriodicBoundaryConditions=true;

long double random;
std::cout << std::setprecision(10) << std::fixed;
long double T,q22,density2;

xim=0.5;
Lsp_13=(xim*L_30)/(1-xim);

etam = 0.7;
Lrg_02=(etam*L_01)/(1-etam);    

xim = 0.6;
Lsp_23=(xim*L_30)/(1-xim);


// q22=1; 
// density2=0.1;
// InitialConditionNonHomogeneous(density2,q22);
// T=200
// SpatialModelSimplified(T);
// cout<<time_sim<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl; 
           
#endif
}

















long double H()
{

                //Co-occurrence Matrix
                long double f11,f12,f21,f22,f;
                int vi,vj;
                f11=0;
                f12=0;
                f21=0;
                f22=0;
                f=0;
                for(int i=0;i<L;i++)
                {
                  for(int j=0;j<L;j++)
                  {
                      if(S[i][j]==1)
                      {
                         //v1
                         vi=i-1;
                         vj=j;
                         if(vi>=0 && vi<L)
                         {
                            if(S[vi][vj]==1){f11++;}
                            if(S[vi][vj]==2){f12++;}
                            f++;
                         }
                         
                         //v2
                         vi=i+1;
                         vj=j;
                         if(vi>=0 && vi<L)
                         {
                            if(S[vi][vj]==1){f11++;}
                            if(S[vi][vj]==2){f12++;}
                            f++;
                         }
                         
                         //v3
                         vi=i;
                         vj=j-1;
                         if(vj>=0 && vj<L)
                         {
                            if(S[vi][vj]==1){f11++;}
                            if(S[vi][vj]==2){f12++;}
                            f++;
                         }
                         
                         //v4
                         vi=i;
                         vj=j+1;
                         if(vj>=0 && vj<L)
                         {
                            if(S[vi][vj]==1){f11++;}
                            if(S[vi][vj]==2){f12++;}
                            f++;
                         }
                      }
                      if(S[i][j]==2)
                      {
                         //v1
                         vi=i-1;
                         vj=j;
                         if(vi>=0 && vi<L)
                         {
                            if(S[vi][vj]==1){f21++;}
                            if(S[vi][vj]==2){f22++;}
                            f++;
                         }
                         
                         //v2
                         vi=i+1;
                         vj=j;
                         if(vi>=0 && vi<L)
                         {
                            if(S[vi][vj]==1){f21++;}
                            if(S[vi][vj]==2){f22++;}
                            f++;
                         }
                         
                         //v3
                         vi=i;
                         vj=j-1;
                         if(vj>=0 && vj<L)
                         {
                            if(S[vi][vj]==1){f21++;}
                            if(S[vi][vj]==2){f22++;}
                            f++;
                         }
                         
                         //v4
                         vi=i;
                         vj=j+1;
                         if(vj>=0 && vj<L)
                         {
                            if(S[vi][vj]==1){f21++;}
                            if(S[vi][vj]==2){f22++;}
                            f++;
                         }
                      }
                  }
                }
                
                long double cmatrix[4];
                cmatrix[0]=f11/f;
                cmatrix[1]=f12/f;
                cmatrix[2]=f21/f;
                cmatrix[3]=f22/f;
                   
                long double entropy;    
                entropy=0;
                for(int i=0;i<4;i++)
                {
                    if(cmatrix[i]!=0)
                    {
                       entropy=entropy-(cmatrix[i]*log2(cmatrix[i]));
                    }
                }
               //cout<<q22<<" "<<entropy<<endl;

    
    return entropy;
}


//---- Grid index mapping (row-major): k = ix*Lp + jx ----
void GridMap(int ix, int jx, int Lp, int &k)
{
    k = ix * Lp + jx;
}

void ReverseGridMap(int &ix, int &jx, int Lp, int k)
{
    ix = k / Lp;
    jx = k % Lp;
}

//---- (Re)allocate S[][] and A[][] for a given grid size, freeing any previous allocation ----
static int allocatedL_S = -1;
static int allocatedN_A = -1;

void FreeGrids()
{
    if (S != nullptr) {
        for (int i = 0; i < allocatedL_S; i++) { delete [] S[i]; }
        delete [] S;
        S = nullptr;
    }
    if (A != nullptr) {
        for (int k = 0; k < allocatedN_A; k++) { delete [] A[k]; }
        delete [] A;
        A = nullptr;
    }
}

void AllocateGrids(int newL)
{
    if (newL == allocatedL_S && S != nullptr && A != nullptr) {
        L = newL;
        N = L * L;
        system_size = N;
        return;
    }
    FreeGrids();
    L = newL;
    N = L * L;
    system_size = N;
    S = new int *[L];
    for (int i = 0; i < L; i++) { S[i] = new int [L]; }
    A = new long double *[N];
    for (int k = 0; k < N; k++) { A[k] = new long double [r]; }
    allocatedL_S = L;
    allocatedN_A = N;
}

//---- RNG seeding: explicit seed for reproducibility, or a time+counter derived seed otherwise ----
static unsigned long seed_call_counter = 0;

void SeedRNG(int seed)
{
    if (seed >= 0) {
        srand((unsigned) seed);
    } else {
        seed_call_counter++;
        unsigned long t = (unsigned long) std::chrono::high_resolution_clock::now().time_since_epoch().count();
        srand((unsigned) (t ^ (seed_call_counter * 2654435761UL)));
    }
}

//---- Fill the entire grid with a single state (e.g. an all-invader or
// all-post-fire-empty base before layering a pattern on top of it). ----
void FillGrid(int state)
{
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { S[i][j]=state; } }
}

//---- Recompute N0..N5/n0..n5 by scanning the whole grid, rather than
// assuming only states 1/2 (the original single-species landscape's only
// possible values) are present. Needed once landscapes can contain any of
// the model's 6 states (multi-layer construction) or are supplied directly
// from R. ----
void TallyGridToCompartments()
{
    N0=0; N1=0; N2=0; N3=0; N4=0; N5=0;
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            switch (S[i][j])
            {
                case 0: N0++; break;
                case 1: N1++; break;
                case 2: N2++; break;
                case 3: N3++; break;
                case 4: N4++; break;
                case 5: N5++; break;
            }
        }
    }
    n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
    n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;
}

void InitialConditionNonHomogeneous(long double &density2, long double &q22, int state1_param, int state2_param, bool reset_grid)
{


    int state1=state1_param;
    int state2=state2_param;


    vector<int> cluster;
    long double nclusters=0;
    int random_int;
    int ix,jx;



    if (reset_grid)
    {
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L;j++)
            {
                S[i][j]=state1;
            }
        }
    }

    //Seed
    // Robustness/generality fix: the original code picked a uniformly
    // random cell out of ALL N cells and just assumed it was state1
    // (always true when reset_grid just filled the whole grid with
    // state1). When reset_grid is false -- placing a second layer's
    // pattern into whatever background cells a prior layer left behind --
    // that assumption can pick a cell already claimed by another state, so
    // we instead collect the actual state1 cells and seed from among
    // those. A single full-grid scan is trivial next to the O(target
    // density * N) growth loop below, and (unlike rejection sampling)
    // can't hang if state1 is sparse.
    {
        vector<int> bg_candidates;
        bg_candidates.reserve(N);
        for (int i=0;i<L;i++)
        {
            for (int j=0;j<L;j++)
            {
                if (S[i][j]==state1) { int idx; GridMap(i,j,L,idx); bg_candidates.push_back(idx); }
            }
        }
        if (bg_candidates.empty())
        {
            Rcpp::stop("InitialConditionNonHomogeneous: no cells in the requested background state are left to seed the pattern from (background fully consumed by a prior layer).");
        }
        random_int = bg_candidates[rand() % bg_candidates.size()];
    }
    ReverseGridMap(ix,jx,L,random_int);
    S[ix][jx]=state2;
    cluster.insert(cluster.end(), random_int);
    nclusters++;
    
    


    
    
    long double threshold1=(state1+state1)*0.5;   //-> 1
    long double threshold2=(state1+state2)*0.5;   //-> 1.5
    long double threshold3=(state2+state2)*0.5;   //-> 2
    
    long double D;
    long double exp11,exp12,exp22;
    long double expr11,expr12,expr22;
    
    
    
    //Calculo de las probabilidades actuales
    long double mean;
    long double est11,est12,est22;
    long double dest11,dest12,dest22;
    int vi,vj;
    est11=0;
    est12=0;
    est22=0;
    
    
    
    //Estimamos p00,p01,p11 inicial
    // Generality fix: when this grid already carries a third state (a
    // prior layer's pattern, or an unrelated background), a neighbor pair
    // touching that state isn't a state1/state2 pair at all and must not
    // be folded into est11/est12/est22 (or its denominator, num_tuples) --
    // otherwise the "similar"/"mixed" classification below, which assumes
    // every pair mean falls at exactly threshold1/2/3, silently miscounts
    // it. In the original single-layer usage the whole grid is always
    // exactly {state1,state2} at this point, so this guard is a no-op
    // there (every pair already qualifies).
    long double num_tuples=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            bool center_in = (S[i][j]==state1 || S[i][j]==state2);

            //Vecino1
            vi=i+1;
            vj=j;
            if(vi>=0 && vi<L && center_in && (S[vi][vj]==state1 || S[vi][vj]==state2))
            {
               mean=(S[vi][vj]+S[i][j])*0.5;
               if(mean<threshold2){est11++;}
               if(mean>threshold1 && mean<threshold3){est12++;}
               if(mean>threshold2){est22++;}
               num_tuples++;
            }

            //Vecino2
            vi=i-1;
            vj=j;
            if(vi>=0 && vi<L && center_in && (S[vi][vj]==state1 || S[vi][vj]==state2))
            {
               mean=(S[vi][vj]+S[i][j])*0.5;
               if(mean<threshold2){est11++;}
               if(mean>threshold1 && mean<threshold3){est12++;}
               if(mean>threshold2){est22++;}
               num_tuples++;
            }

            //Vecino3
            vi=i;
            vj=j+1;
            if(vj>=0 && vj<L && center_in && (S[vi][vj]==state1 || S[vi][vj]==state2))
            {
               mean=(S[vi][vj]+S[i][j])*0.5;
               if(mean<threshold2){est11++;}
               if(mean>threshold1 && mean<threshold3){est12++;}
               if(mean>threshold2){est22++;}
               num_tuples++;
            }

            //Vecino4
            vi=i;
            vj=j-1;
            if(vj>=0 && vj<L && center_in && (S[vi][vj]==state1 || S[vi][vj]==state2))
            {
               mean=(S[vi][vj]+S[i][j])*0.5;
               if(mean<threshold2){est11++;}
               if(mean>threshold1 && mean<threshold3){est12++;}
               if(mean>threshold2){est22++;}
               num_tuples++;
            }

        }
    }
    dest11=est11/num_tuples;
    dest12=est12/num_tuples;
    dest22=est22/num_tuples;
    

    

    int position1,position2;
    int vi_aux,vj_aux;
    int random_int_neigh;
    int zi,zj;
  
  
    long double est11p1,est12p1,est22p1;
    long double est11p2,est12p2,est22p2;  
    long double dest11p1,dest12p1,dest22p1;
    long double dest11p2,dest12p2,dest22p2;
    long double Dp1,Dp2;
    long double mean_ant,mean_pos;

    
    


    // Robustness fix: when this call is one layer of a multi-layer
    // landscape (see generate_landscape_layers_cpp), background_state may
    // hold far fewer cells than density2*N calls for -- e.g. a second
    // layer asking for more cells than a prior layer left available. The
    // relaxed fallbacks below assume *some* background cell always
    // remains to find; once it's truly exhausted, that assumption breaks
    // and the search (and this whole outer loop) would otherwise spin
    // forever. bg_remaining tracks how many state1 cells are actually
    // left (one O(N) scan up front, then decremented per placement below)
    // so the outer loop can stop and warn instead of hanging when the
    // requested density can't be reached.
    long bg_remaining = 0;
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { if (S[i][j]==state1) bg_remaining++; } }
    bool warned_bg_exhausted = false;

    bool relaxed1=false;
    bool relaxed2=false;
    while(cluster.size()<density2*N-1)
    {
       if (bg_remaining <= 0)
       {
           if (!warned_bg_exhausted)
           {
               warned_bg_exhausted = true;
               Rcpp::Rcout << "[ForestFireR] generate_landscape: stopping this layer early -- "
                              "the background state ran out (only " << cluster.size() <<
                              " of the requested " << (long)(density2*N-1) <<
                              " cells could be placed). This happens when a layer's target "
                              "density leaves less background than a later layer needs; "
                              "check your generate_landscape_layers() density budget." << std::endl;
           }
           break;
       }

       //Seleccionamos un vacio y un estado=1 pegado al aglomerado
       //
       // Robustness fix: the original algorithm here is unbounded random
       // rejection sampling for an "isolated" state1 site (no state2
       // neighbor). For several (density2, p) combinations well inside the
       // paper's own stated ranges (e.g. density2=0.5, p=0.5) the pool of
       // such isolated sites empties out long before the target density is
       // reached, and this loop never terminates -- confirmed to hang
       // indefinitely (not just "slow") across grid sizes L=10..100. Since
       // this is a landscape-generation *setup* step, not part of the
       // reaction dynamics itself, we cap the search and fall back to a
       // relaxed selection (drop the "isolated" requirement, or as a last
       // resort accept any remaining state1 site) rather than hang forever.
       // A one-time notice is printed via Rcpp::Rcout when this triggers.
       bool condition1=false;
       long attempts1=0;
       const long MAX_ATTEMPTS1 = 2000L;  // fixed, independent of N -- see note above the outer loop
       while(condition1==false)
       {
        attempts1++;
        if(attempts1 > MAX_ATTEMPTS1)
        {
            if(!relaxed1)
            {
                relaxed1=true;
                Rcpp::Rcout << "[ForestFireR] generate_landscape: relaxing the "
                               "'isolated site' constraint for this replacement "
                               "(no such site found after " << MAX_ATTEMPTS1 <<
                               " attempts) -- this density2/p combination hits an "
                               "algorithmic limitation of the inherited landscape "
                               "generator. Result may deviate slightly from a "
                               "strict Appendix B run." << std::endl;
            }
            // relaxed: accept ANY remaining state1 site, no isolation requirement
            random_int = (rand() % (N));
            ReverseGridMap(ix,jx,L,random_int);
            if(S[ix][jx]==state1)
            {
               condition1=true;
               position1=random_int;
            }
            continue;
        }

        random_int = (rand() % (N));
        ReverseGridMap(ix,jx,L,random_int);
        if(S[ix][jx]==state1)
        {
            int vecinos=1;
            
            vi=ix+1;
            vj=jx;
            if(vi>=0 && vi<L){ if(S[vi][vj]==state2){vecinos=vecinos*0;}}
            vi=ix-1;
            vj=jx;            
            if(vi>=0 && vi<L){ if(S[vi][vj]==state2){vecinos=vecinos*0;}}
            vi=ix;
            vj=jx+1;
            if(vj>=0 && vj<L){ if(S[vi][vj]==state2){vecinos=vecinos*0;}}
            vi=ix;
            vj=jx-1;
            if(vj>=0 && vj<L){ if(S[vi][vj]==state2){vecinos=vecinos*0;}}
            
            if(vecinos==1)
            {
               condition1=true;
               position1=random_int;
            }       
       } 
      }
      bool condition2=false;
      long attempts2=0;
      const long MAX_ATTEMPTS2 = 2000L;  // fixed; its fallback is an O(N) linear scan, so keep this small
      while(condition2==false)
      {
        attempts2++;
        if(attempts2 > MAX_ATTEMPTS2)
        {
            if(!relaxed2)
            {
                relaxed2=true;
                Rcpp::Rcout << "[ForestFireR] generate_landscape: falling back to "
                               "a linear scan for a cluster-adjacent state1 site "
                               "(random search exhausted after " << MAX_ATTEMPTS2 <<
                               " attempts)." << std::endl;
            }
            bool found=false;
            for(int si=0; si<L && !found; si++)
            {
                for(int sj=0; sj<L && !found; sj++)
                {
                    if(S[si][sj]==state1)
                    {
                        int di[4]={1,-1,0,0}; int dj[4]={0,0,1,-1};
                        for(int dd=0; dd<4; dd++)
                        {
                            int ni=si+di[dd], nj=sj+dj[dd];
                            if(ni>=0 && ni<L && nj>=0 && nj<L && S[ni][nj]==state2)
                            {
                                condition2=true;
                                GridMap(si,sj,L,position2);
                                found=true;
                                break;
                            }
                        }
                    }
                }
            }
            if(!found)
            {
                // ultimate fallback: no state1 site borders the cluster at all
                // (only possible if state1 has been fully consumed elsewhere);
                // accept any remaining state1 site so the loop still terminates.
                for(int si=0; si<L && !condition2; si++)
                {
                    for(int sj=0; sj<L && !condition2; sj++)
                    {
                        if(S[si][sj]==state1)
                        {
                            condition2=true;
                            GridMap(si,sj,L,position2);
                        }
                    }
                }
            }
            continue;
        }

          random_int = (rand() % (cluster.size()));
          ReverseGridMap(ix,jx,L,cluster[random_int]);
          random_int_neigh = (rand() % (4));

                       
          if(random_int_neigh==0)
          {
             vi_aux=ix+1;
             vj_aux=jx;
          }
          if(random_int_neigh==1)
          {
             vi_aux=ix-1;
             vj_aux=jx;
          }
          if(random_int_neigh==2)
          {
             vi_aux=ix;
             vj_aux=jx+1;
          }
          if(random_int_neigh==3)
          {
             vi_aux=ix;
             vj_aux=jx-1;
          }
          if(vi_aux>=0 && vi_aux<L && vj_aux>=0 && vj_aux<L)
          {
             if(S[vi_aux][vj_aux]==state1)
             {
                condition2=true;
                GridMap(vi_aux,vj_aux,L,position2);
             }
          }          
      }       
        

      
    long double density0=(cluster.size()+1)/system_size;
    //Expected from Desired Distribution and the actual number of nodes in the cluster
    exp11=1+density0*(q22-2);
    exp12=2*density0*(1-q22);
    exp22=density0*q22;
      
      
     //---------- //Hacemos los cambios de las probabilidades para (vi,vj) lejos del cluster
     // (Both this block and the "cerca cluster" block below carry the same
     // state1/state2 membership guard added to the num_tuples loop above,
     // for the same reason: a neighbor belonging to some third state --
     // another layer's already-placed pattern -- isn't part of this
     // layer's state1/state2 pair statistics and must be skipped, not
     // treated as if it were state1.)
           ReverseGridMap(vi,vj,L,position1);
           est11p1=est11;
           est12p1=est12;
           est22p1=est22;
       
           zi=vi+1;
           zj=vj;
           if(zi>=0 && zi<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p1--;}
             if(mean_pos<threshold2){est11p1++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p1--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p1++;}
             if(mean_ant>threshold2){est22p1--;}
             if(mean_pos>threshold2){est22p1++;}
           }  
           
           zi=vi-1;
           zj=vj;
           if(zi>=0 && zi<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p1--;}
             if(mean_pos<threshold2){est11p1++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p1--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p1++;}
             if(mean_ant>threshold2){est22p1--;}
             if(mean_pos>threshold2){est22p1++;}
           }  
 
           zi=vi;
           zj=vj+1;
           if(zj>=0 && zj<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p1--;}
             if(mean_pos<threshold2){est11p1++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p1--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p1++;}
             if(mean_ant>threshold2){est22p1--;}
             if(mean_pos>threshold2){est22p1++;}
           }  
                     
           zi=vi;
           zj=vj-1;
           if(zj>=0 && zj<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p1--;}
             if(mean_pos<threshold2){est11p1++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p1--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p1++;}
             if(mean_ant>threshold2){est22p1--;}
             if(mean_pos>threshold2){est22p1++;}
           }                      
           dest11p1=est11p1/num_tuples;
           dest12p1=est12p1/num_tuples;
           dest22p1=est22p1/num_tuples;
      
      
           Dp1=fabs(dest11p1-exp11)+fabs(dest12p1-exp12)+fabs(dest22p1-exp22);
      

      
     //---------- //Hacmos los cambios de las probabilidades para (vi,vj) cerca cluster
           ReverseGridMap(vi,vj,L,position2);  
           est11p2=est11;
           est12p2=est12;
           est22p2=est22;
       
           zi=vi+1;
           zj=vj;
           if(zi>=0 && zi<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p2--;}
             if(mean_pos<threshold2){est11p2++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p2--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p2++;}
             if(mean_ant>threshold2){est22p2--;}
             if(mean_pos>threshold2){est22p2++;}
           }  
           
           zi=vi-1;
           zj=vj;
           if(zi>=0 && zi<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p2--;}
             if(mean_pos<threshold2){est11p2++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p2--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p2++;}
             if(mean_ant>threshold2){est22p2--;}
             if(mean_pos>threshold2){est22p2++;}
           }  
 
           zi=vi;
           zj=vj+1;
           if(zj>=0 && zj<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p2--;}
             if(mean_pos<threshold2){est11p2++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p2--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p2++;}
             if(mean_ant>threshold2){est22p2--;}
             if(mean_pos>threshold2){est22p2++;}
           }  
                     
           zi=vi;
           zj=vj-1;
           if(zj>=0 && zj<L && (S[zi][zj]==state1 || S[zi][zj]==state2))
           {
              mean_ant=(S[zi][zj]+state1)*0.5;
              mean_pos=(S[zi][zj]+state2)*0.5;
                            
             //Variacion de las probabilidades ante el nuevo cambio de S[i][j]
             if(mean_ant<threshold2){est11p2--;}
             if(mean_pos<threshold2){est11p2++;}
             if(mean_ant>threshold1 && mean_ant<threshold3){est12p2--;}
             if(mean_pos>threshold1 && mean_pos<threshold3){est12p2++;}
             if(mean_ant>threshold2){est22p2--;}
             if(mean_pos>threshold2){est22p2++;}
           }                      
           dest11p2=est11p2/num_tuples;
           dest12p2=est12p2/num_tuples;
           dest22p2=est22p2/num_tuples;
      
      
           Dp2=fabs(dest11p2-exp11)+fabs(dest12p2-exp12)+fabs(dest22p2-exp22);
      

           // Robustness fix: this was originally two independent `if`s
           // (Dp1<Dp2 / Dp2<Dp1), so an exact tie (Dp1==Dp2) placed
           // nothing and the outer loop made zero progress that
           // iteration. Ties are rare in the "far"/"near" candidates the
           // strict algorithm finds, but became common once the
           // isolated-site/cluster-adjacent-site searches above started
           // falling back to relaxed selection (position1 and position2
           // can then coincide) -- confirmed to cause the outer growth
           // loop to spin forever at high target densities (e.g.
           // density2=0.9) regardless of grid size. Using if/else
           // guarantees exactly one placement happens every iteration.
           if(Dp1<Dp2)
           {
             ReverseGridMap(vi,vj,L,position1);
             S[vi][vj]=state2;
             cluster.insert(cluster.end(),position1);
             bg_remaining--;

             est11=est11p1;
             est12=est12p1;
             est22=est22p1;
           }
           else
           {
             ReverseGridMap(vi,vj,L,position2);
             S[vi][vj]=state2;
             cluster.insert(cluster.end(),position2);
             bg_remaining--;

             est11=est11p2;
             est12=est12p2;
             est22=est22p2;
           }

    }
     
     
    
// Generality fix: the original code only ever tallied states 1/2 here
// (N0/N3/N4/N5 always stayed 0), which was exactly right as long as the
// grid could only ever contain those two states -- true for the original
// single-layer landscape, no longer true once a landscape can be built
// from several layered calls (e.g. two species over an empty background,
// or fire over an invader background -- see generate_landscape_layers_cpp)
// that leave states 0, 3, 4 and/or 5 on the grid too. TallyGridToCompartments()
// scans and counts all 6 states; for the original single-layer case this
// gives identical N1/N2 (and N0=N3=N4=N5=0, since those states never
// appear there either) -- purely additive, not a behavior change for
// existing callers.
TallyGridToCompartments();

}


















void InitialConditionFromRandomness(long double &density2, int &number_clusters)
{


    int random_int;
    int ix,jx;

N0=0;
N1=0;
N2=0;
N3=0;
N4=0;
N5=0;
for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        S[i][j]=1;
        N1++;
    }
}



    vector<int> cluster; 

    //Seeds
    long double nclusters=0;
    while(nclusters<number_clusters)
    {
        random_int = (rand() % (N));
        int ix,jx;
        ReverseGridMap(ix,jx,L,random_int);
        if(S[ix][jx]==1)
        {
           S[ix][jx]=2;
           N1--;
           N2++;
           nclusters++;
           cluster.insert(cluster.end(), random_int); 
        }
    }
    

    

    
    int kx;
    int random_int_neigh;
    while(cluster.size()<density2*N-1)
    {
    

  
          random_int = (rand() % (cluster.size()));
          ReverseGridMap(ix,jx,L,cluster[random_int]);          
          random_int_neigh = (rand() % (4));
         
        
          if(random_int_neigh==0 && (ix+1)>=0 && (ix+1)<L)
          {
           if(S[ix+1][jx]!=2)
           {
             S[ix+1][jx]=2;
             N1--;
             N2++;
             GridMap(ix+1,jx,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
          if(random_int_neigh==1 && (ix-1)>=0 && (ix-1)<L)
          {
           if(S[ix-1][jx]!=2)
           {
             S[ix-1][jx]=2;
             N1--;
             N2++;
             GridMap(ix-1,jx,L,kx);
             cluster.insert(cluster.end(), kx); 
            }
          }
          if(random_int_neigh==2 && (jx+1)>=0 && (jx+1)<L)
          {
           if(S[ix][jx+1]!=2)
           {
             S[ix][jx+1]=2;
             N1--;
             N2++;
             GridMap(ix,jx+1,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
          if(random_int_neigh==3 && (jx-1)>=0 && (jx-1)<L)
          {
           if(S[ix][jx-1]!=2)
           {
             S[ix][jx-1]=2;
             N1--;
             N2++;
             GridMap(ix,jx-1,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
    
    
    }
    


n0=(N4+N5)/system_size;
n1=N1/system_size;
n2=N2/system_size;
n3=N3/system_size;
n4=N4/system_size;
n5=N5/system_size;


}



































long double Entropy(long double &density2)
{

  long double gscale = N*density2;
  int ngrid = sqrt(gscale); 
  int nxgrid=L-ngrid+1;  
  
  int npoints=ngrid*ngrid;  //number of points inside each grid cell
  int ncells=nxgrid*nxgrid; //total number of cells
  
  
  long double sumprob=0;
  long double frequencies[ncells];
  int contador=0;
  for(int i=0;i<nxgrid;i++)
  {
      for(int j=0;j<nxgrid;j++)
      {
          //We count probability in this subsquare
          long double freq=0;
          for(int ki=i;ki<i+ngrid;ki++)
          {
              for(int kj=j;kj<j+ngrid;kj++)
              {
                  if(S[ki][kj]==2){freq++;}
              }
          }
          freq=freq/npoints;
          frequencies[contador]=freq;
          sumprob=sumprob+freq;
          contador++;
      }
  }
  
  long double entropy=0;
  
  for(int i=0;i<ncells;i++)
  {
      frequencies[i]=frequencies[i]/sumprob;
      if(frequencies[i]!=0){entropy=entropy-frequencies[i]*log(frequencies[i]);}
  }
  
  entropy=entropy/log(ncells);
  
  return entropy;
  //cout<<entropy<<endl;

}



















void InitialCondition(long double &density2)
{

     

N0=0;
N1=0;
N2=0;
N3=0;
N4=0;
N5=0;
for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        S[i][j]=1;
        N1++;
    }
}
    
int random_int;
int ix,jx;
for(int i=1;i<density2*N;i++)
{
    random_int = (rand() % (N));
    ReverseGridMap(ix,jx,L,random_int);
    S[ix][jx]=2;
    N1--;
    N2++;
}



n0=(N4+N5)/system_size;
n1=N1/system_size;
n2=N2/system_size;
n3=N3/system_size;
n4=N4/system_size;
n5=N5/system_size;


}





void SpatialModelSimplified(long double &T)
{


long double time,tau;  
long double sumA,sumt;
int ix,jx;
long double random;
int candidatek;
int candidatex,candidatey;
int nx,ny,nk;
int dt_movie=100;
long double paso_anterior=0;
bool with_fire;

std::cout << std::setprecision(10) << std::fixed;

    int pasos=0;

    sumA=0;
    for(int k=0;k<N;k++)
    {     
        VegetationFireModel_Rates(k); 
        for(int j=0;j<r;j++){sumA=sumA+A[k][j];}                 
    }
    if(sumA<0.00000001){sumA=0;}       

    
    bool ext=false;
    int vuelta=0;
    time_sim=0;
    int contador_mv=0;
    int contt=0;
    while(time_sim<T && sumA!=0 && ext==false)
    {

    
    
        n0=(N4+N5)/system_size;
        n1=N1/system_size;
        n2=N2/system_size;
        n3=N3/system_size;
        n4=N4/system_size;
        n5=N5/system_size;
        

        if(n1<0.0001){ext=true;}
        //cout<<time_sim<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl;  

    
        
        random = ((double) rand() / (RAND_MAX));
        tau = -log(random)/sumA;
        time_sim=time_sim+tau;
        
        

        
        random = ((double) rand() / (RAND_MAX));
        sumt=0;  
        bool update=false;
        int k=0;
        while(update==false && sumA!=0)
        {

            for(int j=0;j<r;j++)
            {
                if(sumt/sumA<=random && random<(sumt+A[k][j])/sumA)
                {

                   ReverseGridMap(ix,jx,L,k);
                   if(j==0)
                   {
                         N4--;
                         N1++;
                      
                      S[ix][jx]=1;
                   }           
                   if(j==1)
                   {
                      S[ix][jx]=4;
                      N1--;
                      N4++;
                   }
                   if(j==2)
                   {
                         N4--;
                         N2++;
                      S[ix][jx]=2;
                   }
                   if(j==3)
                   {
                      S[ix][jx]=4;
                      N2--;
                      N4++;
                   }  
                   if(j==4)
                   {
                      S[ix][jx]=2;
                      N1--;
                      N2++;
                   }   
                   if(j==5)
                   {
                      S[ix][jx]=1;
                      N2--;
                      N1++;
                   } 
                   if(j==6)
                   {
                      S[ix][jx]=3;
                      N1--;
                      N3++;
                   }
                   if(j==7)
                   {
                      S[ix][jx]=3;
                      N2--;
                      N3++;
                   }  
                   if(j==8)
                   {
                      S[ix][jx]=5;
                      N3--;
                      N5++;                     
                   }  
                   if(j==9)
                   {
                      S[ix][jx]=1;
                      N5--;
                      N1++;                     
                   } 
                   if(j==10)
                   {
                      S[ix][jx]=2;
                      N5--;
                      N2++;                     
                   } 
                   if(j==11)
                   {
                         N5--;
                         N1++;
                      S[ix][jx]=1;
                   }  
                   if(j==12)
                   {
                         N5--;
                         N2++;
                      S[ix][jx]=2;
                   }  
                   candidatex=ix;
                   candidatey=jx;
                   candidatek=k;
                   
                   update=true; 
          
                }
                sumt=sumt+A[k][j];
            }             
        k++;
        }  
        
       
        
        
    //Actualize Rates of (i,j) and its neighbourhood    
    long double sum_aux=sumA;
    //n0
    nx=candidatex;
    ny=candidatey;
    nk=candidatek;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {        
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n1
    nx=candidatex+1;
    ny=candidatey;
    nk=candidatek+L;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n2
    nx=candidatex-1;
    ny=candidatey;
    nk=candidatek-L;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {          
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n3
    nx=candidatex;
    ny=candidatey+1;
    nk=candidatek+1;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {        
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n4
    nx=candidatex;
    ny=candidatey-1;
    nk=candidatek-1;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {                                                  
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    
    
    
    
    
    //Periodic Boundary Conditions
    
  if(PeriodicBoundaryConditions==true)
  {
    nx=candidatex+1;
    ny=candidatey;   
    if(nx<0)
    {
       nx=nx+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(nx>L-1)
    {
       nx=nx-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    
    nx=candidatex-1;
    ny=candidatey;   
    if(nx<0)
    {
       nx=nx+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(nx>L-1)
    {
       nx=nx-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    nx=candidatex;
    ny=candidatey+1;   
    if(ny<0)
    {
       ny=ny+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(ny>L-1)
    {
       ny=ny-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    nx=candidatex;
    ny=candidatey-1;   
    if(ny<0)
    {
       ny=ny+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(ny>L-1)
    {
       ny=ny-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    //Fin actualizacion de las tasas asociadas a las Boundary Conditions   
  }



          

    
    
    sumA=sum_aux;
    if(sumA<0.00000001){sumA=0;} 


    vuelta++;
    }  
    

  
}
















bool Boundaries(int &ix, int &jx)
{
     bool bounded=true;
     if(ix<0 || ix>L-1 || jx<0 || jx>L-1){bounded=false;}
     return bounded;
}



void InitialConditionInClusters(long double &density2)
{

//We consider only one cluster with no preferential direction of attachment

N0=0;
N1=0;
N2=0;
N3=0;
N4=0;
N5=0;
for(int i=0;i<L;i++)
{
    for(int j=0;j<L;j++)
    {
        S[i][j]=1;
        N1++;
    }
}



    vector<int> cluster; 

    //Seed of V2
    int random_int = (rand() % (N));
    int ix,jx;
    ReverseGridMap(ix,jx,L,random_int);
    S[ix][jx]=2;
    N1--;
    N2++;
    cluster.insert(cluster.end(), random_int); 
    

    
    int kx;
    int random_int_neigh;
    while(cluster.size()<density2*N-1)
    {
    

  
          random_int = (rand() % (cluster.size()));
          ReverseGridMap(ix,jx,L,cluster[random_int]);          
          random_int_neigh = (rand() % (4));
         
        
          if(random_int_neigh==0 && (ix+1)>=0 && (ix+1)<L)
          {
           if(S[ix+1][jx]!=2)
           {
             S[ix+1][jx]=2;
             N1--;
             N2++;
             GridMap(ix+1,jx,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
          if(random_int_neigh==1 && (ix-1)>=0 && (ix-1)<L)
          {
           if(S[ix-1][jx]!=2)
           {
             S[ix-1][jx]=2;
             N1--;
             N2++;
             GridMap(ix-1,jx,L,kx);
             cluster.insert(cluster.end(), kx); 
            }
          }
          if(random_int_neigh==2 && (jx+1)>=0 && (jx+1)<L)
          {
           if(S[ix][jx+1]!=2)
           {
             S[ix][jx+1]=2;
             N1--;
             N2++;
             GridMap(ix,jx+1,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
          if(random_int_neigh==3 && (jx-1)>=0 && (jx-1)<L)
          {
           if(S[ix][jx-1]!=2)
           {
             S[ix][jx-1]=2;
             N1--;
             N2++;
             GridMap(ix,jx-1,L,kx);
             cluster.insert(cluster.end(), kx); 
           }
          }
    
    
    }
    


n0=(N4+N5)/system_size;
n1=N1/system_size;
n2=N2/system_size;
n3=N3/system_size;
n4=N4/system_size;
n5=N5/system_size;


}























void SpatialModel_Movie(long double &T)
{


long double time,tau;  
long double sumA,sumt;
int ix,jx;
long double random;
int candidatek;
int candidatex,candidatey;
int nx,ny,nk;


long double dt_movie=0.1;


   
   
   

    sumA=0;
    for(int k=0;k<N;k++)
    {     
        VegetationFireModel_Rates(k); 
        for(int j=0;j<r;j++){sumA=sumA+A[k][j];}                 
    }
    if(sumA<0.00000001){sumA=0;}       

    
    long double paso_anterior=0;
    int vuelta=0;
    int contador_mv=0;
    time=0;
    do{

    
    
        n0=(N4+N5)/system_size;
        n1=N1/system_size;
        n2=N2/system_size;
        n3=N3/system_size;
        n4=N4/system_size;
        n5=N5/system_size;
     
     



        cout<<time<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl;  
        
        random = ((double) rand() / (RAND_MAX));
        tau = -log(random)/sumA;
        time=time+tau;

        
        random = ((double) rand() / (RAND_MAX));
        sumt=0;  
        bool update=false;
        int k=0;
        do{

            for(int j=0;j<r;j++)
            {
                if(sumt/sumA<=random && random<(sumt+A[k][j])/sumA)
                {

                   ReverseGridMap(ix,jx,L,k);
                   if(j==0)
                   {
                         N4--;
                         N1++;
                      
                      S[ix][jx]=1;
                   }           
                   if(j==1)
                   {
                      S[ix][jx]=4;
                      N1--;
                      N4++;
                   }
                   if(j==2)
                   {
                         N4--;
                         N2++;
                      S[ix][jx]=2;
                   }
                   if(j==3)
                   {
                      S[ix][jx]=4;
                      N2--;
                      N4++;
                   }  
                   if(j==4)
                   {
                      S[ix][jx]=2;
                      N1--;
                      N2++;
                   }   
                   if(j==5)
                   {
                      S[ix][jx]=1;
                      N2--;
                      N1++;
                   } 
                   if(j==6)
                   {
                      S[ix][jx]=3;
                      N1--;
                      N3++;
                   }
                   if(j==7)
                   {
                      S[ix][jx]=3;
                      N2--;
                      N3++;
                   }  
                   if(j==8)
                   {
                      S[ix][jx]=5;
                      N3--;
                      N5++;                     
                   }  
                   if(j==9)
                   {
                      S[ix][jx]=1;
                      N5--;
                      N1++;                     
                   } 
                   if(j==10)
                   {
                      S[ix][jx]=2;
                      N5--;
                      N2++;                     
                   } 
                   if(j==11)
                   {
                         N5--;
                         N1++;
                      S[ix][jx]=1;
                   }  
                   if(j==12)
                   {
                         N5--;
                         N2++;
                      S[ix][jx]=2;
                   }  
                   candidatex=ix;
                   candidatey=jx;
                   candidatek=k;
                   
                   update=true; 
          
                }
                sumt=sumt+A[k][j];
            }             
        k++;
        }while(update==false && sumA!=0);  
        
       
        
        
    //Actualize Rates of (i,j) and its neighbourhood    
    long double sum_aux=sumA;
    //n0
    nx=candidatex;
    ny=candidatey;
    nk=candidatek;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {        
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n1
    nx=candidatex+1;
    ny=candidatey;
    nk=candidatek+L;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n2
    nx=candidatex-1;
    ny=candidatey;
    nk=candidatek-L;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {          
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n3
    nx=candidatex;
    ny=candidatey+1;
    nk=candidatek+1;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {        
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    //n4
    nx=candidatex;
    ny=candidatey-1;
    nk=candidatek-1;
    if(nx>=0 && nx<L && ny>=0 && ny<L)
    {                                                  
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}       
    }
    
    
    
    
    
    //Periodic Boundary Conditions
    
  if(PeriodicBoundaryConditions==true)
  {
    nx=candidatex+1;
    ny=candidatey;   
    if(nx<0)
    {
       nx=nx+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(nx>L-1)
    {
       nx=nx-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    
    nx=candidatex-1;
    ny=candidatey;   
    if(nx<0)
    {
       nx=nx+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(nx>L-1)
    {
       nx=nx-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    nx=candidatex;
    ny=candidatey+1;   
    if(ny<0)
    {
       ny=ny+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(ny>L-1)
    {
       ny=ny-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    
    nx=candidatex;
    ny=candidatey-1;   
    if(ny<0)
    {
       ny=ny+L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    if(ny>L-1)
    {
       ny=ny-L;
       GridMap(nx,ny,L,nk);  
       
       for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
       VegetationFireModel_Rates(nk);
       for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}   
    }
    //Fin actualizacion de las tasas asociadas a las Boundary Conditions   
  }



          

    
    
    sumA=sum_aux;
    if(sumA<0.00000001){sumA=0;} 


    vuelta++;
    }while(time<T && sumA!=0);   
    
    
    
      

    
}



























vector<int> Cluster(int &state)
{


     
   
     int **plantilla; 
     
     plantilla = new int *[L];
     for (int i = 0; i < L; i++){plantilla[i] = new int [L];}

     int nums=0;
     for(int i=0;i<L;i++)
     {
         for(int j=0;j<L;j++)
         {
             plantilla[i][j]=0;
             if(S[i][j]==state)
             {
                plantilla[i][j]=1; 
                nums++;
             }
         }
     }          
    
    
    
  vector<int> clusters;  
  int num_seeds=nums;
  int ncommunity=2;
     
   while(nums!=0)
   {
     
     vector<int> ncluster;  
          
             
     int si,sj;
     int ix,jx;
     bool seed=false;
     int cont=0;
     while(seed==false)
     {
           ReverseGridMap(ix,jx,L,cont);
           if(plantilla[ix][jx]==1)
           {
              si=ix;
              sj=jx;
              seed=true;
              plantilla[si][sj]=ncommunity;
              ncluster.push_back(cont);  
           }           
     cont++;
     }


     int point;
     int pointb;
     int puntero=0;
     do
     {    
         point=ncluster[puntero];
         ReverseGridMap(ix,jx,L,point); 
         
         //buscamos todos los vecinos de point
         if(ix+1>=0 && ix+1<L)
         {
             if(plantilla[ix+1][jx]==1)
             {
                GridMap(ix+1,jx,L,pointb);
                ncluster.push_back(pointb); 
                plantilla[ix+1][jx]=ncommunity;
                
             }
         }
         if(ix+1>=0 && ix+1<L && jx-1>=0 && jx-1<L)
         {
            if(plantilla[ix+1][jx-1]==1)
            {
               GridMap(ix+1,jx-1,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix+1][jx-1]=ncommunity;
            }
         }
         if(ix+1>=0 && ix+1<L && jx+1>=0 && jx+1<L)
         {
            if(plantilla[ix+1][jx+1]==1)
            {
               GridMap(ix+1,jx+1,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix+1][jx+1]=ncommunity;
            }
         }
         if(ix-1>=0 && ix-1<L)
         {
            if(plantilla[ix-1][jx]==1)
            {
               GridMap(ix-1,jx,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix-1][jx]=ncommunity;
            }
         }
         if(ix-1>=0 && ix-1<L && jx-1>=0 && jx-1<L)
         {
            if(plantilla[ix-1][jx-1]==1)
            {
               GridMap(ix-1,jx-1,L,pointb);
               ncluster.push_back(pointb); 
               plantilla[ix-1][jx-1]=ncommunity; 
            }
         }
         if(ix-1>=0 && ix-1<L && jx+1>=0 && jx+1<L)
         {
            if(plantilla[ix-1][jx+1]==1)
            {
               GridMap(ix-1,jx+1,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix-1][jx+1]=ncommunity; 
            }
         }
         if(jx+1>=0 && jx+1<L)
         {         
            if(plantilla[ix][jx+1]==1)
            {
               GridMap(ix,jx+1,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix][jx+1]=ncommunity;
            } 
         } 
         if(jx-1>=0 && jx-1<L)
         {      
            if(plantilla[ix][jx-1]==1)
            {
               GridMap(ix,jx-1,L,pointb);
               ncluster.push_back(pointb);  
               plantilla[ix][jx-1]=ncommunity;
            } 
         } 
                                     
    puntero++;
    }while(puntero<ncluster.size());
           
    clusters.push_back(ncluster.size());  
    ncommunity++;

     
  nums=nums-ncluster.size();
  }
     
     

        
     return clusters;
       
}








void VegetationFireModel_Rates(int &k)  
{
     //This function returns four Reaction Rates (j=1,...,4) for each grid node k. Total number of states=3. Total number of reactions=4.
     //Input: Node Number: k
     //Output: Reaction Rates: A[k][j], for j=1,...,4
     //Notes: It takes into account boundary conditions.
     //Actualiza las tasas A[k][j] del nodo asociado a (ix,jx)
     
        
     long double v1,v2,v3;
     int s;
     int ix,jx;
     int state;
    
     for(int j=0;j<r;j++){A[k][j]=0;}
     ReverseGridMap(ix,jx,L,k);
     
     
          if(PeriodicBoundaryConditions==true)
          {     
            state=1;
            v1=neighbours_PeriodicBoundaryConditions(ix,jx,state);  
               
            state=2;
            v2=neighbours_PeriodicBoundaryConditions(ix,jx,state);
               
            state=3;
            v3=neighbours_PeriodicBoundaryConditions(ix,jx,state);
         }
         
          if(PeriodicBoundaryConditions==false)
          {     
            state=1;
            v1=neighbours(ix,jx,state);  
               
            state=2;
            v2=neighbours(ix,jx,state);
               
            state=3;
            v3=neighbours(ix,jx,state);
         }
            
               
            if(S[ix][jx]==4)
            {                                             
               A[k][0]=L_01*(v1)+Lr_01;                           
               A[k][2]=L_02*(v2)+Lr_02;
            }
            if(S[ix][jx]==1)
            {        
               A[k][1] = L_10;              
               A[k][4] = L_12*(v2)+Lr_12;                                               
               A[k][6] = Lig_13 + Lsp_13*(v3);              
            }
            if(S[ix][jx]==2)
            {        
               A[k][3] = L_20;
               A[k][5] = L_21*(v1)+Lr_21;                                 
               A[k][7] = Lig_23 + Lsp_23*(v3);              
            }
            if(S[ix][jx]==3)
            {
               A[k][8]=L_30;
            }            
            if(S[ix][jx]==5)
            {
               A[k][9]=Lrg_01;
               A[k][10]=Lrg_02;
               A[k][11]=L_01*(v1)+Lr_01;                           
               A[k][12]=L_02*(v2)+Lr_02;
            }
            

}













void MeanField(long double &T)
{
     Vector particle;
     particle.x=n0;  //0, Empty    
     particle.y=n1;  //1, V1
     particle.z=n2;  //2, V2
     particle.q=n3;  //3, F   
     particle.s=n4;  //3, empty_v
     particle.l=n5;  //3, empty_f
   
     particle.t=0;
     //for(double k=0;k<T;k=k+h)

     cout<<particle.t<<" "<<particle.y<<" "<<particle.z<<" "<<particle.q<<" "<<particle.x<<endl;  //Full trajectory        
     long double pasot=0.1;
     do
     {
         particle=zRK4(particle);         
         if(particle.t>=pasot)
         {
            cout<<particle.t<<" "<<particle.y<<" "<<particle.z<<" "<<particle.q<<" "<<particle.x<<endl;  //Full trajectory
            pasot=pasot+0.1;
         }

     }while(particle.t<T);
     //cout<<parameter<<" "<<particle.t<<" "<<particle.y<<" "<<particle.z<<" "<<particle.q<<" "<<particle.x<<endl;  //Full trajectory 
     
     
}






Vector zRK4(Vector &particle)  
{

   //Parameters numerical resolution   
   Vector p;
   Vector k1,k2,k3,k4;

 
   
   //k1
   p.t=particle.t;
   //p.x=particle.x;
   p.y=particle.y;
   p.z=particle.z;
   p.q=particle.q;
   p.s=particle.s;
   p.l=particle.l;
   k1=ForestFire(p);
   

   long double h=0.01;

   
   
   //k2
   p.t=particle.t+0.5*h;
   //p.x=particle.x+0.5*k1.x*h;
   p.y=particle.y+0.5*k1.y*h;
   p.z=particle.z+0.5*k1.z*h;
   p.q=particle.q+0.5*k1.q*h;
   p.s=particle.s+0.5*k1.s*h;
   p.l=particle.l+0.5*k1.l*h;
   k2=ForestFire(p);
   
   //k3
   p.t=particle.t+0.5*h;
   //p.x=particle.x+0.5*k2.x*h;
   p.y=particle.y+0.5*k2.y*h;
   p.z=particle.z+0.5*k2.z*h;
   p.q=particle.q+0.5*k2.q*h;
   p.s=particle.s+0.5*k2.s*h;
   p.l=particle.l+0.5*k2.l*h;
   k3=ForestFire(p);
   
   //K4
   p=particle;
   p.t=particle.t+h;
   //p.x=particle.x+k3.x*h;
   p.y=particle.y+k3.y*h;
   p.z=particle.z+k3.z*h;
   p.q=particle.q+k3.q*h;
   p.s=particle.s+k3.s*h;
   p.l=particle.l+k3.l*h;
   k4=ForestFire(p);
     

   //particle.x=particle.x+h*((k1.x+2*k2.x+2*k3.x+k4.x)/6);
   particle.y=particle.y+h*((k1.y+2*k2.y+2*k3.y+k4.y)/6);
   particle.z=particle.z+h*((k1.z+2*k2.z+2*k3.z+k4.z)/6);
   particle.q=particle.q+h*((k1.q+2*k2.q+2*k3.q+k4.q)/6);
   particle.s=particle.s+h*((k1.s+2*k2.s+2*k3.s+k4.s)/6);
   particle.l=particle.l+h*((k1.l+2*k2.l+2*k3.l+k4.l)/6);
   particle.t=particle.t+h;
   
   particle.x=particle.s+particle.l;
   
   
   
   
      
   return particle;
   

}

























Vector ForestFire(Vector &particle)  
{  
   long double x,y,z,q,t,s,l;
   x=particle.x;   
   y=particle.y;
   z=particle.z;
   q=particle.q;
   s=particle.s;
   l=particle.l;
   t=particle.t;
   
   Vector f;
   
   
     
     
//x=n0 //0, empty   
//y=n1 //V1
//z=n2 //V2
//q=n3 //F
//s=n4 //Empty_v
//l=n5 //Empty_f
     
     
     

     
   f.y = -L_10*y - (L_12*z)*y + (L_21*y)*z  + L_01*y*(s+l) + Lr_01*(s+l) - (Lig_13+Lsp_13*q)*y + (Lrg_01)*l; //V1
   f.z = -L_20*z - (L_21*y)*z + (L_12*z)*y  + L_02*z*(s+l) + Lr_02*(s+l) - (Lig_23+Lsp_23*q)*z + (Lrg_02)*l; //V2   
   f.q = -L_30*q + (Lig_13+Lsp_13*q)*y + (Lig_23+Lsp_23*q)*z; //F
   f.s = (L_10*y + L_20*z) - (L_01*y+L_02*z)*s - (Lr_01 + Lr_02)*s; //Empty_V
   f.l = - (L_01*y+L_02*z)*l - (Lr_01 + Lr_02)*l - (Lrg_01+Lrg_02)*l + q*(L_30); //Empty_F
   
     
   return f;
}








long double neighbours_PeriodicBoundaryConditions(int &ix,int &jx,int &state)
{
     long double v=0;
            

     int ix0=ix-1;
     if(ix0<0){ix0=ix0+L;}
     if(ix0>L-1){ix0=ix0-L;}
     
     int ix1=ix+1;
     if(ix1<0){ix1=ix1+L;}
     if(ix1>L-1){ix1=ix1-L;}
     
     int jx0=jx-1;
     if(jx0<0){jx0=jx0+L;}
     if(jx0>L-1){jx0=jx0-L;}
     
     int jx1=jx+1;
     if(jx1<0){jx1=jx1+L;}
     if(jx1>L-1){jx1=jx1-L;}
     
     
     
     if(S[ix0][jx]==state){v++;}
     if(S[ix1][jx]==state){v++;}
     if(S[ix][jx0]==state){v++;}
     if(S[ix][jx1]==state){v++;}  
     
       

       return v; 
}











void MeanField_Guillespie(long double &T)
{


long double random;
long double time,tau;
int kr=r;
long double W[kr];
long double sumA=1;


   
    time=0;
    do{

        n0=(N4+N5)/system_size;
        n1=N1/system_size;
        n2=N2/system_size;
        n3=N3/system_size;
        n4=N4/system_size;
        n5=N5/system_size;
        
                         
        W[0]=N4*(L_01*N1);
        W[1]=N1*(L_10);
        W[2]=N4*(L_02*N2);    
        W[3]=N2*(L_20);
        W[4]=N1*(L_12*N2);
        W[5]=N2*(L_21*N1);
        W[6]=N1*(Lsp_13*N3+Lig_13);
        W[7]=N2*(Lsp_23*N3+Lig_23);
        W[8]=N3*(L_30);
        W[9]=N5*(Lrg_01);
        W[10]=N5*(Lrg_02);        
        W[11]=N5*(L_01*N1);
        W[12]=N5*(L_02*N2);
        
        
        sumA=0;
        for(int i=0;i<kr;i++){sumA=sumA+W[i];}

        random = ((double) rand() / (RAND_MAX));
        tau = -log(random)/sumA;
        time=time+tau;
        cout<<time<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl;                       


        random = ((double) rand() / (RAND_MAX));           
        if(random<(W[0]/sumA))
        {          
           N4--;
           N1++;
        }    
        if(random>=(W[0]/sumA) && random<((W[0]+W[1])/sumA))
        {
           N1--;
           N4++; 
        }
        if(random>=((W[0]+W[1])/sumA) && random<((W[0]+W[1]+W[2])/sumA))
        {
           N4--;
           N2++;
        }
        if(random>=((W[0]+W[1]+W[2])/sumA) && random<((W[0]+W[1]+W[2]+W[3])/sumA))
        {
           N2--;
           N4++; 
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4])/sumA))
        {
           N1--;
           N2++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5])/sumA))
        {
           N2--;
           N1++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6])/sumA))
        {
           N1--;
           N3++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7])/sumA))
        {
           N2--;
           N3++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8])/sumA))
        {
           N3--;
           N5++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9])/sumA))
        {
           N5--;
           N1++;           
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9]+W[10])/sumA))
        {
           N5--;
           N2++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9]+W[10])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9]+W[10]+W[11])/sumA))
        {
           N5--;
           N1++;
        } 
        if(random>=((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9]+W[10]+W[11])/sumA) && random<((W[0]+W[1]+W[2]+W[3]+W[4]+W[5]+W[6]+W[7]+W[8]+W[9]+W[10]+W[11]+W[12])/sumA))
        {
           N5--;
           N2++;
        } 
                         
    }while(time<T && sumA!=0);


}





long double neighbours(int &ix,int &jx,int &state)
{
     long double v=0;
            
     if(ix>0 && ix<L-1 && jx>0 && jx<L-1) 
     {    
            if(S[ix-1][jx]==state){v++;}
            if(S[ix+1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;}
            if(S[ix][jx+1]==state){v++;}  
     }
       
     //Boundary Conditions
     if(ix==0)
     {
        if(jx==0)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx+1]==state){v++;}                     
        }        
        if(jx==L-1)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx-1]==state){v++;}
        }     
        if(jx>0 && jx<L-1)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx-1]==state){v++;}
           if(S[ix][jx+1]==state){v++;}
        }          
      } 
      if(ix==L-1)
      {
         if(jx==0)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx+1]==state){v++;}
         }                
         if(jx==L-1)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;}                 
         }
         if(jx>0 && jx<L-1)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;}
            if(S[ix][jx+1]==state){v++;}              
         }
       } 
       if(jx==0 && ix>0 && ix<L-1)
       {
          if(S[ix-1][jx]==state){v++;}
          if(S[ix+1][jx]==state){v++;}
          if(S[ix][jx+1]==state){v++;}             
       }
       if(jx==L-1 && ix>0 && ix<L-1)
       {
          if(S[ix-1][jx]==state){v++;}
          if(S[ix+1][jx]==state){v++;}
          if(S[ix][jx-1]==state){v++;}
       }    
       return v; 
}




long double full_neighbours(int &ix,int &jx,int &state)
{
     long double v=0;
            
     if(ix>0 && ix<L-1 && jx>0 && jx<L-1) 
     {    
            if(S[ix-1][jx]==state){v++;}
            if(S[ix+1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;}
            if(S[ix][jx+1]==state){v++;}             
            if(S[ix-1][jx-1]==state){v++;}
            if(S[ix-1][jx+1]==state){v++;}
            if(S[ix+1][jx-1]==state){v++;}
            if(S[ix+1][jx+1]==state){v++;}  
     }
       
     //Boundary Conditions
     if(ix==0)
     {
        if(jx==0)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx+1]==state){v++;} 
           if(S[ix+1][jx+1]==state){v++;}                     
        }        
        if(jx==L-1)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx-1]==state){v++;}
           if(S[ix+1][jx-1]==state){v++;}                     
        }     
        if(jx>0 && jx<L-1)
        {
           if(S[ix+1][jx]==state){v++;}
           if(S[ix][jx-1]==state){v++;}
           if(S[ix][jx+1]==state){v++;}
           if(S[ix+1][jx-1]==state){v++;}                     
           if(S[ix+1][jx+1]==state){v++;}                     
        }          
      } 
      if(ix==L-1)
      {
         if(jx==0)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx+1]==state){v++;}
            if(S[ix-1][jx+1]==state){v++;}
         }                
         if(jx==L-1)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;} 
            if(S[ix-1][jx-1]==state){v++;}                
         }
         if(jx>0 && jx<L-1)
         {
            if(S[ix-1][jx]==state){v++;}
            if(S[ix][jx-1]==state){v++;}
            if(S[ix][jx+1]==state){v++;}      
            if(S[ix-1][jx+1]==state){v++;}
            if(S[ix-1][jx-1]==state){v++;}        
         }
       } 
       if(jx==0 && ix>0 && ix<L-1)
       {
          if(S[ix-1][jx]==state){v++;}
          if(S[ix+1][jx]==state){v++;}
          if(S[ix][jx+1]==state){v++;} 
          if(S[ix-1][jx+1]==state){v++;}  
          if(S[ix+1][jx+1]==state){v++;}              
       }
       if(jx==L-1 && ix>0 && ix<L-1)
       {
          if(S[ix-1][jx]==state){v++;}
          if(S[ix+1][jx]==state){v++;}
          if(S[ix][jx-1]==state){v++;}
          if(S[ix-1][jx-1]==state){v++;}
          if(S[ix+1][jx-1]==state){v++;}
       }    
       return v;
}








void SpatialModel(long double &T)
{

long double time,tau;  
long double sumA,sumt;
int ix,jx;
long double random;



    time=0;
    while(time<T)
    {
    
        sumA=0;
        for(int k=0;k<N;k++)
        {     
            VegetationFireModel_Rates(k); 
            for(int j=0;j<r;j++)
            {
                sumA=sumA+A[k][j];
            }
        } 
        
        

        n0=(N4+N5)/system_size;
        n1=N1/system_size;
        n2=N2/system_size;
        n3=N3/system_size;
        n4=N4/system_size;
        n5=N5/system_size;
        cout<<time<<" "<<n1<<" "<<n2<<" "<<n3<<" "<<n0<<endl;  
        
        

        
        
        
        
        random = ((double) rand() / (RAND_MAX));
        tau = -log(random)/sumA;
        
      
        random = ((double) rand() / (RAND_MAX));
        sumt=0;
        for(int k=0;k<N;k++)
        {
            for(int j=0;j<r;j++)
            {
                if(sumt<=random*sumA && random*sumA<(sumt+A[k][j]))
                {   

                   //transicion en k y transicion j
                   ReverseGridMap(ix,jx,L,k);
                   if(j==0)
                   {
                      S[ix][jx]=1;
                      N4--;
                      N1++;                     
                   }           
                   if(j==1)
                   {
                      S[ix][jx]=4;
                      N1--;
                      N4++;
                   }
                   if(j==2)
                   {
                         N4--;
                         N2++;
                      S[ix][jx]=2;
                   }
                   if(j==3)
                   {
                      S[ix][jx]=4;
                      N2--;
                      N4++;
                   }  
                   if(j==4)
                   {
                      S[ix][jx]=2;
                      N1--;
                      N2++;
                   }   
                   if(j==5)
                   {
                      S[ix][jx]=1;
                      N2--;
                      N1++;
                   } 
                   if(j==6)
                   {
                      S[ix][jx]=3;
                      N1--;
                      N3++;
                   }
                   if(j==7)
                   {
                      S[ix][jx]=3;
                      N2--;
                      N3++;
                   }  
                   if(j==8)
                   {
                      S[ix][jx]=5;
                      N3--;
                      N5++;                     
                   }  
                   if(j==9)
                   {
                      S[ix][jx]=1;
                      N5--;
                      N1++;                     
                   } 
                   if(j==10)
                   {
                      S[ix][jx]=2;
                      N5--;
                      N2++;                     
                   } 
                   if(j==11)
                   {
                         N5--;
                         N1++;
                      S[ix][jx]=1;
                   }  
                   if(j==12)
                   {
                         N5--;
                         N2++;
                      S[ix][jx]=2;
                   }   
                } 
                
            sumt=sumt+A[k][j];              
            }
        }

      
        
    time=time+tau;       
    }

}


//=====================================================================
// Trajectory-recording globals (shared by the *_Rec functions below)
//=====================================================================
vector<long double> rec_time, rec_n0, rec_n1, rec_n2, rec_n3, rec_n4, rec_n5;

static void ClearRecording()
{
    rec_time.clear();
    rec_n0.clear(); rec_n1.clear(); rec_n2.clear();
    rec_n3.clear(); rec_n4.clear(); rec_n5.clear();
}

static void PushRecording(long double t)
{
    rec_time.push_back(t);
    rec_n0.push_back(n0); rec_n1.push_back(n1); rec_n2.push_back(n2);
    rec_n3.push_back(n3); rec_n4.push_back(n4); rec_n5.push_back(n5);
}

//=====================================================================
// SpatialModelSimplified_Rec: same Gillespie SSA as SpatialModelSimplified,
// but optionally records the (time, n0..n5) trajectory instead of only
// keeping the final densities.
//   record_dt < 0  -> no trajectory recorded (fastest, only final state)
//   record_dt == 0 -> record every accepted event (dense, memory heavy)
//   record_dt > 0  -> record a sample every record_dt time units
//=====================================================================
void SpatialModelSimplified_Rec(long double &T, long double &record_dt, bool check_extinction)
{
    ClearRecording();

    long double tau;
    long double sumA, sumt;
    int ix, jx;
    long double random;
    int candidatek;
    int candidatex, candidatey;
    int nx, ny, nk;
    long double next_record = 0;

    sumA = 0;
    for (int k = 0; k < N; k++) {
        VegetationFireModel_Rates(k);
        for (int j = 0; j < r; j++) { sumA = sumA + A[k][j]; }
    }
    if (sumA < 0.00000001) { sumA = 0; }

    bool ext = false;
    time_sim = 0;

    n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
    n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;
    if (record_dt >= 0) { PushRecording(time_sim); next_record = record_dt; }

    while (time_sim < T && sumA != 0 && ext == false) {

        n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
        n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;

        if (check_extinction && n1 < 0.0001) { ext = true; }

        random = ((double) rand() / (RAND_MAX));
        tau = -log(random) / sumA;
        time_sim = time_sim + tau;

        random = ((double) rand() / (RAND_MAX));
        sumt = 0;
        bool update = false;
        int k = 0;
        while (update == false && sumA != 0) {
            for (int j = 0; j < r; j++) {
                if (sumt / sumA <= random && random < (sumt + A[k][j]) / sumA) {
                    ReverseGridMap(ix, jx, L, k);
                    if (j==0)  { N4--; N1++; S[ix][jx]=1; }
                    if (j==1)  { S[ix][jx]=4; N1--; N4++; }
                    if (j==2)  { N4--; N2++; S[ix][jx]=2; }
                    if (j==3)  { S[ix][jx]=4; N2--; N4++; }
                    if (j==4)  { S[ix][jx]=2; N1--; N2++; }
                    if (j==5)  { S[ix][jx]=1; N2--; N1++; }
                    if (j==6)  { S[ix][jx]=3; N1--; N3++; }
                    if (j==7)  { S[ix][jx]=3; N2--; N3++; }
                    if (j==8)  { S[ix][jx]=5; N3--; N5++; }
                    if (j==9)  { S[ix][jx]=1; N5--; N1++; }
                    if (j==10) { S[ix][jx]=2; N5--; N2++; }
                    if (j==11) { N5--; N1++; S[ix][jx]=1; }
                    if (j==12) { N5--; N2++; S[ix][jx]=2; }
                    candidatex = ix; candidatey = jx; candidatek = k;
                    update = true;
                }
                sumt = sumt + A[k][j];
            }
            k++;
        }

        long double sum_aux = sumA;
        nx=candidatex; ny=candidatey; nk=candidatek;
        if (nx>=0 && nx<L && ny>=0 && ny<L) {
            for (int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
            VegetationFireModel_Rates(nk);
            for (int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}
        }
        nx=candidatex+1; ny=candidatey; nk=candidatek+L;
        if (nx>=0 && nx<L && ny>=0 && ny<L) {
            for (int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
            VegetationFireModel_Rates(nk);
            for (int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}
        }
        nx=candidatex-1; ny=candidatey; nk=candidatek-L;
        if (nx>=0 && nx<L && ny>=0 && ny<L) {
            for (int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
            VegetationFireModel_Rates(nk);
            for (int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}
        }
        nx=candidatex; ny=candidatey+1; nk=candidatek+1;
        if (nx>=0 && nx<L && ny>=0 && ny<L) {
            for (int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
            VegetationFireModel_Rates(nk);
            for (int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}
        }
        nx=candidatex; ny=candidatey-1; nk=candidatek-1;
        if (nx>=0 && nx<L && ny>=0 && ny<L) {
            for (int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];}
            VegetationFireModel_Rates(nk);
            for (int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];}
        }

        if (PeriodicBoundaryConditions == true) {
            nx=candidatex+1; ny=candidatey;
            if (nx<0)   { nx=nx+L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }
            if (nx>L-1) { nx=nx-L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }

            nx=candidatex-1; ny=candidatey;
            if (nx<0)   { nx=nx+L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }
            if (nx>L-1) { nx=nx-L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }

            nx=candidatex; ny=candidatey+1;
            if (ny<0)   { ny=ny+L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }
            if (ny>L-1) { ny=ny-L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }

            nx=candidatex; ny=candidatey-1;
            if (ny<0)   { ny=ny+L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }
            if (ny>L-1) { ny=ny-L; GridMap(nx,ny,L,nk); for(int j=0;j<r;j++){sum_aux=sum_aux-A[nk][j];} VegetationFireModel_Rates(nk); for(int j=0;j<r;j++){sum_aux=sum_aux+A[nk][j];} }
        }

        sumA = sum_aux;
        if (sumA < 0.00000001) { sumA = 0; }

        n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
        n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;

        if (record_dt > 0) {
            while (time_sim >= next_record) { PushRecording(time_sim); next_record += record_dt; }
        } else if (record_dt == 0) {
            PushRecording(time_sim);
        }
    }

    n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
    n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;
    if (record_dt >= 0) { PushRecording(time_sim); }
}

//=====================================================================
// MeanField_Rec: deterministic RK4 integration of the mean-field ODEs
// (Eq. 12 in the paper), recording the full trajectory at intervals of dt.
//=====================================================================
//---- Slow-subsystem RHS: same equations as ForestFire() for y,z,s,l, but with
//     the fire compartment q held FIXED (it is integrated separately, see below). ----
struct SlowState { long double y, z, s, l; };

SlowState ForestFireSlow(long double y, long double z, long double q, long double s, long double l)
{
    SlowState f;
    f.y = -L_10*y - (L_12*z)*y + (L_21*y)*z  + L_01*y*(s+l) + Lr_01*(s+l) - (Lig_13+Lsp_13*q)*y + (Lrg_01)*l;
    f.z = -L_20*z - (L_21*y)*z + (L_12*z)*y  + L_02*z*(s+l) + Lr_02*(s+l) - (Lig_23+Lsp_23*q)*z + (Lrg_02)*l;
    f.s = (L_10*y + L_20*z) - (L_01*y+L_02*z)*s - (Lr_01 + Lr_02)*s;
    f.l = - (L_01*y+L_02*z)*l - (Lr_01 + Lr_02)*l - (Lrg_01+Lrg_02)*l + q*(L_30);
    return f;
}

void MeanField_Rec(long double &T, long double &dt, double n1_0, double n2_0, double n3_0, double n4_0, double n5_0)
{
    // Integrates the same mean-field equations as MeanField()/zRK4()/ForestFire(),
    // but with an IMEX (implicit-explicit) splitting instead of plain explicit RK4:
    // the fire compartment q obeys dq/dt = -k*q + c with k = L_30 - Lsp_13*y - Lsp_23*z,
    // c = Lig_13*y + Lig_23*z. With the paper's default L_30 ~ 1e6 (fire burns out almost
    // instantly) this term is extremely stiff and a fixed-step explicit RK4 on q alone
    // diverges to NaN for any h large enough to be computationally practical. Since the
    // ODE for q is linear given (y,z) frozen over a short step, it has a closed-form
    // solution; using that exact update for q (unconditionally stable for any h) and then
    // integrating the comparatively slow y,z,s,l block explicitly with classic RK4 removes
    // the instability without changing the model equations themselves.
    ClearRecording();

    long double y=n1_0, z=n2_0, q=n3_0, s=n4_0, l=n5_0, x=s+l, t=0;
    long double h = 0.01; // internal step for the slow (y,z,s,l) block

    rec_time.push_back(t);
    rec_n0.push_back(x); rec_n1.push_back(y); rec_n2.push_back(z);
    rec_n3.push_back(q); rec_n4.push_back(s); rec_n5.push_back(l);

    long double pasot = dt;
    do {
        // 1) exact update of the fast fire compartment q, holding y,z fixed over this step
        long double k = L_30 - Lsp_13*y - Lsp_23*z;
        long double c = Lig_13*y + Lig_23*z;
        long double q_new;
        if (k > 1e-8L) {
            q_new = q*expl(-k*h) + (c/k)*(1.0L - expl(-k*h));
        } else if (k < -1e-8L) {
            // k<0 means the linearized fire term is locally growing (spread
            // outrunning extinction) -- can only happen transiently/at
            // extreme parameter corners outside the paper's calibrated
            // range; guard against blow-up with a bounded explicit step.
            q_new = q + (c - k*q)*h;
        } else {
            q_new = q + c*h;
        }
        if (!std::isfinite((double)q_new) || q_new < 0) { q_new = 0; }
        if (q_new > 1) { q_new = 1; }

        // 2) classic RK4 for the slow block (y,z,s,l), holding q at its just-updated value
        SlowState k1 = ForestFireSlow(y, z, q_new, s, l);
        SlowState k2 = ForestFireSlow(y+0.5*h*k1.y, z+0.5*h*k1.z, q_new, s+0.5*h*k1.s, l+0.5*h*k1.l);
        SlowState k3 = ForestFireSlow(y+0.5*h*k2.y, z+0.5*h*k2.z, q_new, s+0.5*h*k2.s, l+0.5*h*k2.l);
        SlowState k4 = ForestFireSlow(y+h*k3.y, z+h*k3.z, q_new, s+h*k3.s, l+h*k3.l);

        y = y + h*((k1.y+2*k2.y+2*k3.y+k4.y)/6);
        z = z + h*((k1.z+2*k2.z+2*k3.z+k4.z)/6);
        s = s + h*((k1.s+2*k2.s+2*k3.s+k4.s)/6);
        l = l + h*((k1.l+2*k2.l+2*k3.l+k4.l)/6);
        q = q_new;
        t = t + h;
        x = s + l;

        if (!std::isfinite((double)y) || y<0){y=0;} if (y>1){y=1;}
        if (!std::isfinite((double)z) || z<0){z=0;} if (z>1){z=1;}
        if (!std::isfinite((double)s) || s<0){s=0;} if (s>1){s=1;}
        if (!std::isfinite((double)l) || l<0){l=0;} if (l>1){l=1;}

        if (t >= pasot) {
            rec_time.push_back(t);
            rec_n0.push_back(x); rec_n1.push_back(y); rec_n2.push_back(z);
            rec_n3.push_back(q); rec_n4.push_back(s); rec_n5.push_back(l);
            pasot = pasot + dt;
        }
    } while (t < T);

    n0=x; n1=y; n2=z; n3=q; n4=s; n5=l;
    time_sim = t;
}

//=====================================================================
// MeanFieldGillespie_Rec: well-mixed (non-spatial) stochastic Gillespie
// simulation, recording the trajectory. Expects N1..N5 and system_size
// to already be set by the caller (see simulate_mean_field_stochastic_cpp).
//=====================================================================
void MeanFieldGillespie_Rec(long double &T, long double &record_dt, long double &sysN)
{
    ClearRecording();
    system_size = sysN;

    long double random, time, tau;
    int kr = r;
    long double W[13];
    long double sumA;
    long double next_record = 0;

    time = 0;
    n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
    n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;
    if (record_dt >= 0) { PushRecording(time); next_record = record_dt; }

    do {
        n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
        n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;

        W[0]=N4*(L_01*N1);
        W[1]=N1*(L_10);
        W[2]=N4*(L_02*N2);
        W[3]=N2*(L_20);
        W[4]=N1*(L_12*N2);
        W[5]=N2*(L_21*N1);
        W[6]=N1*(Lsp_13*N3+Lig_13);
        W[7]=N2*(Lsp_23*N3+Lig_23);
        W[8]=N3*(L_30);
        W[9]=N5*(Lrg_01);
        W[10]=N5*(Lrg_02);
        W[11]=N5*(L_01*N1);
        W[12]=N5*(L_02*N2);

        sumA=0;
        for (int i=0;i<kr;i++){sumA=sumA+W[i];}
        if (sumA < 0.00000001) { sumA = 0; break; }

        random = ((double) rand() / (RAND_MAX));
        tau = -log(random)/sumA;
        time = time + tau;

        random = ((double) rand() / (RAND_MAX));
        long double cum = 0;
        int chosen = -1;
        for (int i=0;i<kr;i++) {
            cum += W[i];
            if (random < cum/sumA) { chosen = i; break; }
        }
        if (chosen < 0) { chosen = kr-1; }

        switch(chosen) {
            case 0:  N4--; N1++; break;
            case 1:  N1--; N4++; break;
            case 2:  N4--; N2++; break;
            case 3:  N2--; N4++; break;
            case 4:  N1--; N2++; break;
            case 5:  N2--; N1++; break;
            case 6:  N1--; N3++; break;
            case 7:  N2--; N3++; break;
            case 8:  N3--; N5++; break;
            case 9:  N5--; N1++; break;
            case 10: N5--; N2++; break;
            case 11: N5--; N1++; break;
            case 12: N5--; N2++; break;
        }

        if (record_dt > 0) {
            while (time >= next_record) { PushRecording(time); next_record += record_dt; }
        } else if (record_dt == 0) {
            PushRecording(time);
        }

    } while (time < T && sumA != 0);

    n0=(N4+N5)/system_size; n1=N1/system_size; n2=N2/system_size;
    n3=N3/system_size; n4=N4/system_size; n5=N5/system_size;
    time_sim = time;
    if (record_dt >= 0) { PushRecording(time); }
}

//=====================================================================
// Rcpp export layer
//=====================================================================

// [[Rcpp::export]]
void set_seed_cpp(int seed) {
    SeedRNG(seed);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix generate_landscape_cpp(int Lgrid, double density2, double p, int seed) {
    AllocateGrids(Lgrid);
    SeedRNG(seed);
    long double density2_ld = density2;
    long double p_ld = p;
    InitialConditionNonHomogeneous(density2_ld, p_ld);
    Rcpp::IntegerMatrix grid(L, L);
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { grid(i,j) = S[i][j]; } }
    return grid;
}

// [[Rcpp::export]]
Rcpp::List simulate_spatial_cpp(
    double T, int Lgrid, double density2, double p,
    double L_01_, double L_02_, double L_10_, double L_20_,
    double L_12_, double L_21_,
    double L_30_, double Lig_13_, double Lig_23_,
    double Lsp_13_, double Lsp_23_,
    double Lrg_01_, double Lrg_02_,
    double Lr_01_, double Lr_02_, double Lr_12_, double Lr_21_,
    bool periodic, int seed,
    double record_dt, bool record_grid)
{
    L_01=L_01_; L_02=L_02_; L_10=L_10_; L_20=L_20_;
    L_12=L_12_; L_21=L_21_; L_30=L_30_;
    Lig_13=Lig_13_; Lig_23=Lig_23_;
    Lsp_13=Lsp_13_; Lsp_23=Lsp_23_;
    Lrg_01=Lrg_01_; Lrg_02=Lrg_02_;
    Lr_01=Lr_01_; Lr_02=Lr_02_; Lr_12=Lr_12_; Lr_21=Lr_21_;
    PeriodicBoundaryConditions = periodic;

    AllocateGrids(Lgrid);
    SeedRNG(seed);

    long double density2_ld = density2;
    long double p_ld = p;
    InitialConditionNonHomogeneous(density2_ld, p_ld);

    Rcpp::IntegerMatrix initial_grid;
    if (record_grid) {
        initial_grid = Rcpp::IntegerMatrix(L, L);
        for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { initial_grid(i,j) = S[i][j]; } }
    }

    long double T_ld = T;
    long double record_dt_ld = record_dt;
    SpatialModelSimplified_Rec(T_ld, record_dt_ld);

    Rcpp::List result;
    result["time_sim"] = (double) time_sim;
    result["n0"] = (double) n0;
    result["n1"] = (double) n1;
    result["n2"] = (double) n2;
    result["n3"] = (double) n3;
    result["n4"] = (double) n4;
    result["n5"] = (double) n5;

    if (record_dt >= 0) {
        result["trajectory"] = Rcpp::DataFrame::create(
            Rcpp::Named("time") = Rcpp::wrap(rec_time),
            Rcpp::Named("n0") = Rcpp::wrap(rec_n0),
            Rcpp::Named("n1") = Rcpp::wrap(rec_n1),
            Rcpp::Named("n2") = Rcpp::wrap(rec_n2),
            Rcpp::Named("n3") = Rcpp::wrap(rec_n3),
            Rcpp::Named("n4") = Rcpp::wrap(rec_n4),
            Rcpp::Named("n5") = Rcpp::wrap(rec_n5)
        );
    }
    if (record_grid) {
        result["initial_grid"] = initial_grid;
        Rcpp::IntegerMatrix final_grid(L, L);
        for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { final_grid(i,j) = S[i][j]; } }
        result["final_grid"] = final_grid;
    }
    return result;
}

// Builds a landscape from a SEQUENCE of layers instead of a single
// background/pattern pair -- needed for initial conditions the original
// generator can't express, e.g. two species placed independently over a
// shared empty background (Figs. 8-9: 5% native + 5% invader over an
// all-post-fire-empty domain), or a small pattern placed over an
// already-uniform domain (Figs. 10-11: 1% fire over an all-invader
// domain). Grid is first filled entirely with fill_state, then each layer
// k grows pattern_states[k] to densities[k]*L*L cells (as a fraction of
// the WHOLE grid, matching how density2 already works elsewhere) with
// heterogeneity ps[k], drawing candidates only from cells still in
// background_states[k] -- so layers are consumed in order and each one
// only ever overwrites background left by the layers before it. All four
// vector arguments must have the same length (one entry per layer).
// [[Rcpp::export]]
Rcpp::IntegerMatrix generate_landscape_layers_cpp(int Lgrid, int fill_state,
    Rcpp::IntegerVector background_states, Rcpp::IntegerVector pattern_states,
    Rcpp::NumericVector densities, Rcpp::NumericVector ps, int seed)
{
    int nlayers = pattern_states.size();
    if (background_states.size() != nlayers || densities.size() != nlayers || ps.size() != nlayers) {
        Rcpp::stop("generate_landscape_layers_cpp: background_states, pattern_states, densities and ps must all have the same length.");
    }

    AllocateGrids(Lgrid);
    SeedRNG(seed);
    FillGrid(fill_state);

    for (int k=0; k<nlayers; k++) {
        long double density_ld = densities[k];
        long double p_ld = ps[k];
        InitialConditionNonHomogeneous(density_ld, p_ld, background_states[k], pattern_states[k], false);
    }

    Rcpp::IntegerMatrix grid(L, L);
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { grid(i,j) = S[i][j]; } }
    return grid;
}

// Runs the spatial Gillespie dynamics from a landscape supplied directly by
// R (typically from generate_landscape_layers_cpp) instead of building one
// internally from density2/p -- the counterpart simulate_spatial_cpp needs
// for any initial condition beyond its own single background/pattern
// generator, e.g. the Fig. 8-11 comparisons.
// [[Rcpp::export]]
Rcpp::List simulate_spatial_from_grid_cpp(
    double T, Rcpp::IntegerMatrix initial_grid,
    double L_01_, double L_02_, double L_10_, double L_20_,
    double L_12_, double L_21_,
    double L_30_, double Lig_13_, double Lig_23_,
    double Lsp_13_, double Lsp_23_,
    double Lrg_01_, double Lrg_02_,
    double Lr_01_, double Lr_02_, double Lr_12_, double Lr_21_,
    bool periodic, int seed,
    double record_dt, bool record_grid, bool check_extinction = true)
{
    L_01=L_01_; L_02=L_02_; L_10=L_10_; L_20=L_20_;
    L_12=L_12_; L_21=L_21_; L_30=L_30_;
    Lig_13=Lig_13_; Lig_23=Lig_23_;
    Lsp_13=Lsp_13_; Lsp_23=Lsp_23_;
    Lrg_01=Lrg_01_; Lrg_02=Lrg_02_;
    Lr_01=Lr_01_; Lr_02=Lr_02_; Lr_12=Lr_12_; Lr_21=Lr_21_;
    PeriodicBoundaryConditions = periodic;

    int Lgrid = initial_grid.nrow();
    if (initial_grid.ncol() != Lgrid) {
        Rcpp::stop("simulate_spatial_from_grid_cpp: initial_grid must be square.");
    }

    AllocateGrids(Lgrid);
    SeedRNG(seed);
    for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { S[i][j] = initial_grid(i,j); } }
    TallyGridToCompartments();

    Rcpp::IntegerMatrix initial_grid_out;
    if (record_grid) {
        initial_grid_out = Rcpp::IntegerMatrix(L, L);
        for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { initial_grid_out(i,j) = S[i][j]; } }
    }

    long double T_ld = T;
    long double record_dt_ld = record_dt;
    SpatialModelSimplified_Rec(T_ld, record_dt_ld, check_extinction);

    Rcpp::List result;
    result["time_sim"] = (double) time_sim;
    result["n0"] = (double) n0;
    result["n1"] = (double) n1;
    result["n2"] = (double) n2;
    result["n3"] = (double) n3;
    result["n4"] = (double) n4;
    result["n5"] = (double) n5;

    if (record_dt >= 0) {
        result["trajectory"] = Rcpp::DataFrame::create(
            Rcpp::Named("time") = Rcpp::wrap(rec_time),
            Rcpp::Named("n0") = Rcpp::wrap(rec_n0),
            Rcpp::Named("n1") = Rcpp::wrap(rec_n1),
            Rcpp::Named("n2") = Rcpp::wrap(rec_n2),
            Rcpp::Named("n3") = Rcpp::wrap(rec_n3),
            Rcpp::Named("n4") = Rcpp::wrap(rec_n4),
            Rcpp::Named("n5") = Rcpp::wrap(rec_n5)
        );
    }
    if (record_grid) {
        result["initial_grid"] = initial_grid_out;
        Rcpp::IntegerMatrix final_grid(L, L);
        for (int i=0;i<L;i++) { for (int j=0;j<L;j++) { final_grid(i,j) = S[i][j]; } }
        result["final_grid"] = final_grid;
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::DataFrame simulate_mean_field_cpp(
    double T, double n1_0, double n2_0, double n3_0, double n4_0, double n5_0,
    double L_01_, double L_02_, double L_10_, double L_20_,
    double L_12_, double L_21_, double L_30_,
    double Lig_13_, double Lig_23_, double Lsp_13_, double Lsp_23_,
    double Lrg_01_, double Lrg_02_,
    double Lr_01_, double Lr_02_, double Lr_12_, double Lr_21_,
    double record_dt)
{
    L_01=L_01_; L_02=L_02_; L_10=L_10_; L_20=L_20_;
    L_12=L_12_; L_21=L_21_; L_30=L_30_;
    Lig_13=Lig_13_; Lig_23=Lig_23_; Lsp_13=Lsp_13_; Lsp_23=Lsp_23_;
    Lrg_01=Lrg_01_; Lrg_02=Lrg_02_;
    Lr_01=Lr_01_; Lr_02=Lr_02_; Lr_12=Lr_12_; Lr_21=Lr_21_;

    long double T_ld = T;
    long double dt_ld = (record_dt > 0) ? record_dt : (T/1000.0);
    MeanField_Rec(T_ld, dt_ld, n1_0, n2_0, n3_0, n4_0, n5_0);

    return Rcpp::DataFrame::create(
        Rcpp::Named("time") = Rcpp::wrap(rec_time),
        Rcpp::Named("n0") = Rcpp::wrap(rec_n0),
        Rcpp::Named("n1") = Rcpp::wrap(rec_n1),
        Rcpp::Named("n2") = Rcpp::wrap(rec_n2),
        Rcpp::Named("n3") = Rcpp::wrap(rec_n3),
        Rcpp::Named("n4") = Rcpp::wrap(rec_n4),
        Rcpp::Named("n5") = Rcpp::wrap(rec_n5)
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame simulate_mean_field_stochastic_cpp(
    double T, double sysN, double n1_0, double n2_0, double n3_0, double n4_0, double n5_0,
    double L_01_, double L_02_, double L_10_, double L_20_, double L_12_, double L_21_, double L_30_,
    double Lig_13_, double Lig_23_, double Lsp_13_, double Lsp_23_, double Lrg_01_, double Lrg_02_,
    int seed, double record_dt)
{
    L_01=L_01_; L_02=L_02_; L_10=L_10_; L_20=L_20_;
    L_12=L_12_; L_21=L_21_; L_30=L_30_;
    Lig_13=Lig_13_; Lig_23=Lig_23_; Lsp_13=Lsp_13_; Lsp_23=Lsp_23_;
    Lrg_01=Lrg_01_; Lrg_02=Lrg_02_;

    SeedRNG(seed);
    system_size = sysN;
    N1 = n1_0 * sysN; N2 = n2_0 * sysN; N3 = n3_0 * sysN;
    N4 = n4_0 * sysN; N5 = n5_0 * sysN;

    long double T_ld = T;
    long double dt_ld = record_dt;
    long double sysN_ld = sysN;
    MeanFieldGillespie_Rec(T_ld, dt_ld, sysN_ld);

    return Rcpp::DataFrame::create(
        Rcpp::Named("time") = Rcpp::wrap(rec_time),
        Rcpp::Named("n0") = Rcpp::wrap(rec_n0),
        Rcpp::Named("n1") = Rcpp::wrap(rec_n1),
        Rcpp::Named("n2") = Rcpp::wrap(rec_n2),
        Rcpp::Named("n3") = Rcpp::wrap(rec_n3),
        Rcpp::Named("n4") = Rcpp::wrap(rec_n4),
        Rcpp::Named("n5") = Rcpp::wrap(rec_n5)
    );
}









