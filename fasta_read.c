#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

#include "fasta.h"

static long
fastaExtraIDfromDescription(const char *fastaIDline)
{
	char *firstBarLocation = NULL;
	long extractedID;

	firstBarLocation = index(fastaIDline, '|');
	if (sscanf(firstBarLocation, "%ld|", &extractedID) != 1) {
		return -1;
	}
	return extractedID;
}

int
fastaReadRecord(FILE *ifp, FASTArecord *fRecord)
{
	/**
	 * A line is "recommended" to be no more than 80 characters,
	 * and the longest sequence is likely no more than 10 lines
	 * of sequence.  The "description" portion is potentially
	 * longer than 80 characters, so we ensure that the relatively
	 * large buffer allocated here is sufficient for that length
	 * as the first action within this function.
	 *
	 * We can allocate a fairly large buffer here without undue
	 * concern that we are wasting memory as this buffer will be
	 * returned to the system as a local variable when this function
	 * returns. */
	char linebuffer[MAX_SEQUENCE_LINES * RECOMMENDED_LINE_LENGTH];
	char *fgetstatus;
	int curLoadIndex = 0, nLinesRead = 0;
	int bytesRemain, curBytesRead;

	/** if our assumption about the first line length is too large
	 * to fit into the allocated buffer, panic
	 */
	assert(MAX_SEQUENCE_LINES * RECOMMENDED_LINE_LENGTH
			> MAX_DESCRIPTION_LINE_LENGTH);

	/** fill linebuffer with zeros */
	bzero(linebuffer, MAX_SEQUENCE_LINES * RECOMMENDED_LINE_LENGTH);

	linebuffer[0] = fgetc(ifp);
	if (linebuffer[0] <= 0) {
		return 0;
	}


	/**
	 * Handle the description.
	 *
	 * Read the rest of the line at the beginning of the loop
	 * adding the result into the current buffer
	 */
	curLoadIndex = 1;
	fgetstatus = fgets(&linebuffer[curLoadIndex],
			MAX_DESCRIPTION_LINE_LENGTH, ifp);
	if (fgetstatus == NULL) {
		fprintf(stderr, "Error: FASTA parser encountered EOF"
				" during partial description line\n");
		printf("RETURNING from %d\n", __LINE__);
		return -1;
	}

	/* check if we have overflow */
	if (linebuffer[strlen(linebuffer) - 1] != '\n') {
		fprintf(stderr, "Error: FASTA parser read description"
				" line greater than %d characters\n",
				MAX_DESCRIPTION_LINE_LENGTH);
		fprintf(stderr, "desc[%s] %lu '%c'\n", linebuffer, strlen(linebuffer),
				linebuffer[strlen(linebuffer)-1]);
		printf("RETURNING from %d\n", __LINE__);
		return -1;
	}
	nLinesRead++;
	fRecord->description = strdup(linebuffer);
	fRecord->id = fastaExtraIDfromDescription(linebuffer);


	/** handle the sequence */
	curLoadIndex = 0;
	bytesRemain = (MAX_SEQUENCE_LINES * RECOMMENDED_LINE_LENGTH) - 2;
	linebuffer[curLoadIndex] = fgetc(ifp);
	if (linebuffer[curLoadIndex] <= 0) {
		fprintf(stderr, "Error: FASTA parser encountered"
					" unexpected end of file before sequence data\n");
		free(fRecord->description);
		printf("RETURNING from %d\n", __LINE__);
		return -1;
	}

	/** collate all of the portions of the sequence */
	while (linebuffer[curLoadIndex] > 0 && linebuffer[curLoadIndex] != '>') {
		fgetstatus = fgets(&linebuffer[curLoadIndex], bytesRemain, ifp);
		if (fgetstatus == NULL) {
			fprintf(stderr, "Error: FASTA parser encountered"
					" unexpected end of file\n");
			free(fRecord->description);
		printf("RETURNING from %d\n", __LINE__);
			return -1;
		}
		nLinesRead++;
		curBytesRead = strlen(&linebuffer[curLoadIndex]);
		if (curBytesRead >= 80) {
			fprintf(stderr,
					"Warning: FASTA parser read sequence of length (%d);",
					curBytesRead);
			fprintf(stderr,
					" longer than 80 character recommendation!\n");
			fprintf(stderr, "       : [%s]\n", &linebuffer[curLoadIndex]);
		}

		if (linebuffer[curLoadIndex + curBytesRead - 1] != '\n') {
			fprintf(stderr, "Error: FASTA parser sequence"
					" line overflows buff\n");
		printf("RETURNING from %d\n", __LINE__);
			return -1;
		}

		/* Now check what the first character of the next line is,
		 * and if it is a new record, push it back */
		bytesRemain -= curBytesRead;
		curLoadIndex += curBytesRead - 1;
		linebuffer[curLoadIndex] = fgetc(ifp);
		if (linebuffer[curLoadIndex] == '>') {
			ungetc('>', ifp);
			linebuffer[curLoadIndex] = 0;
		}
	}

	/** save the sequence */
	fRecord->sequence = strdup(linebuffer);

	return nLinesRead;
}


int
fastaPrintRecord(FILE *ofp, FASTArecord *fRecord)
{
	fprintf(ofp, "FASTA Record:\n");
	fprintf(ofp, "ID   (%ld)\n", fRecord->id);
	fprintf(ofp, "DESC [%s]\n", fRecord->description);
	fprintf(ofp, "SEQ  [%s]\n", fRecord->sequence);

	return 0;
}

/**
 * Allocate and initialize a new FASTA record
 */
void
fastaInitializeRecord(FASTArecord *fRecord)
{
	fRecord->description = NULL;
	fRecord->id = -1;
	fRecord->sequence = NULL;
}

/**
 * Allocate and initialize a new FASTA record
 */
FASTArecord *
fastaAllocateRecord()
{
	FASTArecord *fRecord = NULL;
	
	fRecord = (FASTArecord *) malloc(sizeof(FASTArecord));

	fastaInitializeRecord(fRecord);

	return fRecord;
}

/**
 * clear the record but do not free the fRecord pointer
 */
void
fastaClearRecord(FASTArecord *fRecord)
{
	if (fRecord->description != NULL) {
		free(fRecord->description);
		fRecord->description = NULL;
	}
	if (fRecord->sequence != NULL) {
		free(fRecord->sequence);
		fRecord->sequence = NULL;
	}
	fRecord->id = -1;
}

/**
 * deallocate
 */
void
fastaDeallocateRecord(FASTArecord *fRecord)
{
	fastaClearRecord(fRecord);
	free(fRecord);
}

