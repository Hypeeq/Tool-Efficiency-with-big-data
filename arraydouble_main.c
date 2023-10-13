#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "fasta.h"

// Define an initial size for the array of FASTA records
#define INITIAL_RECORD_ARRAY_SIZE 1

void usage(char *progname)
{
    fprintf(stderr, "%s [<OPTIONS>] <file> [ <file> ...]\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Prints timing of loading and storing FASTA records.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options: \n");
    fprintf(stderr, "-R <REPEATS> : Number of times to repeat load.\n");
    fprintf(stderr, "             : Time reported will be average time.\n");
    fprintf(stderr, "\n");
}

int processFasta(char *filename, FASTArecord **records, double *timeTaken)
{
    FILE *fp;
    FASTArecord fRecord;
    int recordNumber = 0;
    int eofSeen = 0;
    clock_t startTime, endTime;
    int recordArraySize = INITIAL_RECORD_ARRAY_SIZE;
    int currentRecordCount = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failure opening %s: %s\n", filename, strerror(errno));
        return -1;
    }

    // Allocate an initial array for FASTA records
    *records = (FASTArecord *)malloc(recordArraySize * sizeof(FASTArecord));
    if (*records == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        return -1;
    }

    // Record the start time
    startTime = clock();

    do {
        if ((recordNumber % 10000) == 0) {
            printf(".");
            fflush(stdout);
        }

        fastaInitializeRecord(&fRecord);

        int status = fastaReadRecord(fp, &fRecord);
        if (status == 0) {
            eofSeen = 1;
        } else if (status > 0) {
            recordNumber++;
            // Grow the array if necessary
            if (currentRecordCount == recordArraySize) {
                recordArraySize *= 2;
                *records = (FASTArecord *)realloc(*records, recordArraySize * sizeof(FASTArecord));
                if (*records == NULL) {
                    fprintf(stderr, "Memory reallocation error.\n");
                    return -1;
                }
            }
            // Copy the FASTA record to the array
            (*records)[currentRecordCount] = fRecord;
            currentRecordCount++;
            fastaInitializeRecord(&fRecord);
        } else {
            fprintf(stderr, "status = %d\n", status);
            fprintf(stderr, "Error: failure in '%s'\n", filename);
            return -1;
        }

    } while (!eofSeen);
    printf(" %d FASTA records -- %d allocated (%.3lf%% waste)\n", recordNumber, recordArraySize, (((double)(recordArraySize - recordNumber) / recordArraySize) * 100));

    // Record the end time and calculate the difference
    endTime = clock();
    *timeTaken = ((double)(endTime - startTime)) / CLOCKS_PER_SEC;

    fclose(fp);

    return currentRecordCount;
}

int processFastaRepeatedly(char *filename, long repeatsRequested) {
    double timeThisIterationInSeconds;
    double totalTimeInSeconds = 0;
    int minutesPortion;
    int status;
    long i;
    FASTArecord *records = NULL; 

    for (i = 0; i < repeatsRequested; i++) {
        status = processFasta(filename, &records, &timeThisIterationInSeconds);
        if (status < 0) return -1;
        totalTimeInSeconds += timeThisIterationInSeconds;

        // Free memory for FASTA records after processing each iteration
        if (records) {
            for (int j = 0; j < status; j++) {
                free(records[j].description);
                free(records[j].sequence);
            }
            free(records);
            records = NULL;  // Reset the pointer to NULL
        }
    }

    printf("%lf seconds taken for processing total\n", totalTimeInSeconds);

    totalTimeInSeconds /= (double) repeatsRequested;

    minutesPortion = (int) (totalTimeInSeconds / 60);
    totalTimeInSeconds = totalTimeInSeconds - (60 * minutesPortion);
    printf("On average: %d minutes, %lf second per run\n", minutesPortion, totalTimeInSeconds);

    return status;
}


int main(int argc, char **argv)
{
    long repeatsRequested = 1;
    int totalRecordsProcessed = 0;

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            if (argv[i][1] == 'R') {
                if (i + 1 < argc) {
                    if (sscanf(argv[i + 1], "%ld", &repeatsRequested) != 1) {
                        fprintf(stderr, "Error: cannot parse repeats requested from '%s'\n", argv[i + 1]);
                        return 1;
                    }
                    i++;  // Skip the next argument since it has been processed
                } else {
                    fprintf(stderr, "Error: need argument for repeats requested\n");
                    return 1;
                }
            } else {
                fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
                usage(argv[0]);
                return 1;
            }
        } else {
            int recordsProcessed = processFastaRepeatedly(argv[i], repeatsRequested);
            if (recordsProcessed < 0) {
                fprintf(stderr, "Error: Processing '%s' failed -- exiting\n", argv[i]);
                return 1;
            }
            totalRecordsProcessed += recordsProcessed;
        }
    }

    if (totalRecordsProcessed == 0) {
        fprintf(stderr, "No data processed -- provide the name of a file on the command line\n");
        usage(argv[0]);
        return 1;
    }

	

    return 0;
}
