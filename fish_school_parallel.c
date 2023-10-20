// author: Anfernee Pontilan Alviar (22886082), Vicky Chow (23638279)
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>
#include <mpi.h>

#define ROOT_NUM 0;
#define NUM_STEPS 10
#define NUM_FISH 1000
#define FISH_INIT_WEIGHT 15

// Declare structure for fish, holding coordinates (for now)
typedef struct _fish {
    int prev_x;
    int prev_y;
    int x;
    int y;
    double prev_f_i;
    double f_i;
    double delta_f_i;
    double prev_weight;
    double weight;  
} FISH;

double CollectiveAction(FISH* fishes, int num_fish, double total_obj_func) {
    
    double total_distance_times_weight = 0.0;

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < num_fish; i++) 
        {

            double weight = fishes[i].weight;
            double distance = sqrt(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);

            total_distance_times_weight += (distance * weight);

        }
    }

    // Calculate the barycenter
    double barycenter = total_distance_times_weight / total_obj_func;
    return barycenter;
}


// Function for parallelized fish swimming
void swim(FISH* fishes, int num_fish) {
    #pragma omp parallel 
    {
        #pragma omp for
        for (int i = 0; i < num_fish; i++) 
        {   
            int random_x_int = rand() % 201;
            int random_y_int = rand() % 201;
                
            double random_x_movement = (random_x_int / 1000.0) - 0.1;
            double random_y_movement = (random_y_int / 1000.0) - 0.1;
                
            fishes[i].prev_x = fishes[i].x;
            fishes[i].prev_y = fishes[i].y;
                
            fishes[i].x = fishes[i].x + random_x_movement;
            fishes[i].y = fishes[i].y + random_y_movement;
             
        }
    }
}


double get_max_delta(FISH *fishes, int num_fish)
{
    double max_delta = 0.0;

    #pragma omp parallel
    { 
        #pragma omp for // reduction(max:maxDelta)
        for (int i = 0; i < num_fish; i++) 
        {
            fishes[i].delta_f_i = fishes[i].f_i - fishes[i].prev_f_i;

            if (fishes[i].delta_f_i > max_delta) 
            {
                max_delta = fishes[i].delta_f_i;
            }
        }
    }

    return max_delta;
}
// Function for parallelized weight_function
void weight_function(FISH* fishes, int num_fish, double max_delta) 
{
    #pragma omp parallel 
    {   
        #pragma omp for
        for (int i = 0; i < num_fish; i++) {
            fishes[i].weight += (fishes[i].delta_f_i / max_delta);
        }
    }
}

double calc_euc_dist (FISH fish)
{
    return sqrt(pow((double)fish.x, 2) + pow((double)fish.y, 2));
}

double obj_func (FISH* fishes, int num_fish) 
{
    
    double total_sum = 0;

    #pragma omp parallel
    {
        #pragma omp for 
        for (int i = 0; i < num_fish; i++)
        {

            // Set the current f_i to the previous f_i
            fishes[i].prev_f_i = fishes[i].f_i;
            // printf("Fish #%d previous objective function: %f\n", i+1, fishes[i].prev_f_i);

            // Get the distance of the current fish from the center
            double distance = (double)(fishes[i].x * fishes[i].x + fishes[i].y * fishes[i].y);

            // Get the objective function of the current fish
            double obj_func_val = sqrt(distance);
            
            // Store previous objective function value
            fishes[i].prev_f_i = fishes[i].f_i;

            // Set the current objective function for the fish
            fishes[i].f_i = obj_func_val;

            // Add the current f_i to the overall objective function
            total_sum += fishes[i].f_i;
        }
    }   

    return total_sum;
}


int main(int argc, char* argv[])
{
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int rank, num_procs, fishes_per_process;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    FISH *fishes, *local_fishes;

    fishes = (FISH*) malloc(NUM_FISH * sizeof(FISH));
    local_fishes = (FISH*) malloc(fishes_per_process * sizeof(FISH));

    double start = omp_get_wtime();
    
    // Generate positions for the fish; handled by master process
    if (rank == 0)
    {
        #pragma omp parallel for
        for (int i = 0; i < NUM_FISH; i++)
        {
            int x_rand_num = rand() % 201 - 100;
            int y_rand_num = rand() % 201 - 100;

            fishes[i].x = x_rand_num;
            fishes[i].y = y_rand_num;
            fishes[i].weight = FISH_INIT_WEIGHT; 
            fishes[i].f_i = calc_euc_dist(fishes[i]);

        }
    }

    // Determine how many fish are distributed per process
    fishes_per_process = NUM_FISH / num_procs; 

    // Send fish to different processes
    MPI_Scatter(fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 
        local_fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    double local_total_sum, global_total_sum, barycentre;

    for (int i = 0; i < NUM_STEPS; i++)
    {
        
        if (i == 0)
        {
            local_total_sum = obj_func(fishes, fishes_per_process);

            // Perform reduction to calculate global total sum, which will be used to calculate barycentre
            MPI_Reduce(&local_total_sum, &global_total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            
            /*   Perform OpenMP parallelisation on the loop; 
             *   uses local_fishes variable in order to be gathered altogether later 
             */
            #pragma omp parallel 
            {
                #pragma omp for 
                for (int j = 0; j < fishes_per_process; j++)
                {
                    double current_weight = local_fishes[j].weight;
                    // Weight function is random value at the very first step
                    double weight_func = rand() % 5 - 1;
                    local_fishes[j].weight += weight_func;

                }
            }
            
            // Perform swim function on the fishes 
            swim(local_fishes, fishes_per_process);

            // Gather all fish for barycentre calculation
            MPI_Gather(local_fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 
                fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 0, MPI_COMM_WORLD);

            
            // Master process performs collective action
            if(rank == 0)
            {
                barycentre = CollectiveAction(fishes, NUM_FISH, global_total_sum);
            }

            // Information then scattered to worker processes
            MPI_Scatter(fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 
                local_fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 0, MPI_COMM_WORLD);
        }
        else
        {
            // ADD ALL WEIGHTS FOR FISH - LOCAL
            local_total_sum = obj_func(fishes);

            // Get the global total sum across all processes for barycentre calculation
            MPI_Reduce(&local_total_sum, &global_total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            double local_max_delta, global_max_delta;
            
            // Get local max delta for each process
            local_max_delta = get_max_delta(local_fishes, fishes_per_process);

            // Get global max delta to use for calculating new weights
            MPI_Allreduce(&local_max_delta, &global_max_delta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            // Calculate weight for all fishes
            weight_function(local_fishes, fishes_per_process, global_max_delta);  

            // Perform swim function
            swim(local_fishes, fishes_per_process);

            // Gather all fishes
            MPI_Gather(local_fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 
                fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 0, MPI_COMM_WORLD);


            // Perform collective action
            if(rank == 0)
            {
                barycentre = CollectiveAction(fishes, NUM_FISH, global_total_sum);
            }

            // Scatter fish
            MPI_Scatter(fishes, fishes_per_process * sizeof(FISH), MPI_BYTE, 
                local_fishes, fishes_per_process * sizeof(FISH), 0, MPI_COMM_WORLD);
        }
        
    }

    double end = omp_get_wtime();

    double time_taken= end-start;

    free(local_fishes);

    if(rank == 0) 
    {
        free(fishes);
        printf("Time spent: %.2f\n", time_taken);
    }

    MPI_Finalize();
}