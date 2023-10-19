#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h> 
#include <math.h>
#include <stddef.h>

#define NUM_FISH 21
#define FISH_INIT_WEIGHT 15

#define MASTER 0 // task id of first task
#define FROM_MASTER 1 // message type
#define FROM_WORKER 2 // message type

// Declare structure for fish, holding coordinates (for now)
typedef struct _fish {
    int x;
    int y;
} FISH;


int main(int argc, char* argv[])
{
    FISH *fishes;
    int num_tasks; // number of tasks in partition
    int task_id; // task identifier
    int num_workers; // number of worker tasks
    int source; // task id of message source
    int dest; // task id of message destination
    int msg_type; // message type
    int num_fish; // number of fish sent to each worker
    int avg_partition, extra, offset; // determine number of fish sent to each worker
    int i, j, k, rc; // miscellaneous
    MPI_Status status;
    
    // Create type for struct fish
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    const int       nitems=2;
    int             block_len[2] = {1, 1};
    MPI_Datatype    types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype    mpi_fish_type; 
    MPI_Aint        offsets[2];

    offsets[0] = offsetof(FISH, x);
    offsets[1] = offsetof(FISH, y);

    MPI_Type_create_struct(nitems, block_len, offsets, types, &mpi_fish_type);
    MPI_Type_commit(&mpi_fish_type);

    fishes = (FISH*) malloc(NUM_FISH * sizeof(FISH));
    
    if (num_tasks < 2) 
    {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    num_workers = num_tasks;

    if (task_id == MASTER)
    {

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
            fprintf(fp1, "Fish #%d coordinates:  %d, %d\n", i+1, fishes[i].x, fishes[i].y);
        }
        fclose(fp1);

        avg_partition = NUM_FISH / num_workers;
        extra = NUM_FISH % num_workers; 
        offset = 0;
        msg_type = FROM_MASTER; 

        for(dest = 0; dest < num_workers; dest++)
        {
            num_fish = (dest < extra) ? avg_partition+1 : avg_partition;
            printf("Sending %d fish to node %d offset %d\n", num_fish, dest+1, offset);
            MPI_Send(&offset, 1, MPI_INT, dest, msg_type, MPI_COMM_WORLD);
            MPI_Send(&num_fish, 1, MPI_INT, dest, msg_type, MPI_COMM_WORLD);

            MPI_Send(&fishes[offset], num_fish, mpi_fish_type, dest, msg_type, MPI_COMM_WORLD);

            offset += num_fish;
        }

        msg_type = FROM_WORKER;
        for (int i = 0; i < num_workers; i++)
        {
            source = i;

            MPI_Recv(&offset, 1, MPI_INT, source, msg_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&num_fish, 1, MPI_INT, source, msg_type, MPI_COMM_WORLD, &status);

            MPI_Recv(&fishes[offset], num_fish, mpi_fish_type, MASTER, msg_type, MPI_COMM_WORLD, &status);

            printf("Received results from task %d\n", source);

        }

        FILE *fp2;
        fp2 = fopen("post_send.txt", "w+");
        for(int i = 0; i < NUM_FISH; i++)
        {
            fprintf(fp2, "Fish #%d\t", i+1);
            fprintf(fp2, "Coordinates: %d, %d\n", fishes[i].x, fishes[i].y);
        }
        fclose(fp2);

    }

    if (task_id > MASTER)
    {
        msg_type = FROM_MASTER;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, msg_type, MPI_COMM_WORLD, &status);
        MPI_Recv(&num_fish, 1, MPI_INT, MASTER, msg_type, MPI_COMM_WORLD, &status);

        
        MPI_Recv(fishes, num_fish, mpi_fish_type, MASTER, msg_type, MPI_COMM_WORLD, &status);

        for (int i = 0; i < num_fish; i++)
        {
            printf("Fish coordinates: %d, %d\n", fishes[i].x, fishes[i].y);
        }
        
        msg_type = FROM_WORKER;
        MPI_Send(&offset, 1, MPI_INT, MASTER, msg_type, MPI_COMM_WORLD);
        MPI_Send(&num_fish, 1, MPI_INT, MASTER, msg_type, MPI_COMM_WORLD);
       
        MPI_Send(&fishes, num_fish, mpi_fish_type, MASTER, msg_type, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    free(fishes);

    return 0;
}