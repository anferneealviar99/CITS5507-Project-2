#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>

#define NUM_FISH 125
#define FISH_INIT_WEIGHT 15

#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

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


int main(int argc, char* argv[])
{
    FISH *fishes;
    int num_nodes = 4;
    int node_id; // task identifier
    int num_fish; 
    int num_workers;
    int avg_partition; 
    int extras;
    int offset;
    int mtype;
    int source;
    int status;

    fishes = (FISH) malloc(NUM_FISH * sizeof(FISH));

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);

    if (num_nodes < 2)
    {
        printf("Need at least 2 MPI tasks. Aborting...\n")
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(EXIT_FAILURE);
    }
    
    num_workers = num_nodes; 

    if (node_id == MASTER)
    {
        printf("mpi_mm has started with %d nodes.\n", num_nodes);
        printf("Generating fish data...\n");

        for (int i = 0; i < NUM_FISH; i++)
        {
            int x_rand_num = rand() % 201 - 100;
            int y_rand_num = rand() % 201 - 100;

            // Simplified fish init for testing
            fishes[i].x = x_rand_num;
            fishes[i].y = y_rand_num;
        }

        FILE* fp1; 

        fp1 = fopen("pre_send.txt", "w+");
        for(int i = 0; i < NUM_FISH; i++)
        {
            fprintf(fp1, "Fish #%d\t", i+1);
            fprintf(fp1, "Coordinates: %d, %d\n", fishes[i].x, fishes[i].y);
        }

        avg_partition = NUM_FISH / num_workers;
        extras = NUM_FISH % num_workers; 
        offset = 0;
        mtype = FROM_MASTER; 

        for(int dest=0; dest<num_workers; dest++)
        {
            num_fish = (dest <= extras) ? avg_partition+1 : avg_partition;
            printf("Sending %d fish to node %d offset %d\n", num_fish, dest+1, offset)
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            
            for(int i = offset; i < num_fish + offset; i++)
            {
                MPI_Send(&fishes[i].x, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
                MPI_Send(&fishes.[i].y, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            }

            offset = offset + num_fish;
        }

        mtype = FROM_WORKER;
        for (i = 0; i < num_workers; i++)
        {
            source = i;

            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&num_fish, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);

            for (int j = offset; j < num_fish + offset; j++);
            {
                MPI_Recv(&fishes[j].x, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&fishes[j].y, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            }
        }

        FILE *fp2;
        fp2 = fopen("post_send.txt", "w+");
        for(int i = 0; i < NUM_FISH; i++)
        {
            fprintf(fp1, "Fish #%d\t", i+1);
            fprintf(fp1, "Coordinates: %d, %d\n", fishes[i].x, fishes[i].y);
        }

    }

    if (node_id > MASTER)
    {
        mtype = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&num_fish, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);

        for (int j = offset; j < num_fish + offset; j++);
        {
            MPI_Recv(&fishes[j].x, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&fishes[j].y, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
        }
        
        mtype = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&num_fish, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        
        for(int i = offset; i < num_fish + offset; i++)
        {
            MPI_Send(&fishes[i].x, 1, MPI_INT, source, mtype, MPI_COMM_WORLD);
            MPI_Send(&fishes[i].y, 1, MPI_INT, source, mtype, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();

    free(fishes);

    return 0;
}