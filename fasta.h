#ifndef	__FASTA_RECORD_TOOLS_HEADER__
#define	__FASTA_RECORD_TOOLS_HEADER__

typedef struct FASTArecord {
	long id;
	char *description;
	char *sequence;
} FASTArecord;

#define	MAX_SEQUENCE_LINES 1024
#define	RECOMMENDED_LINE_LENGTH 80
#define	MAX_DESCRIPTION_LINE_LENGTH 1024

int  fastaReadRecord(FILE *ifp, FASTArecord *fRecord);
void fastaInitializeRecord(FASTArecord *fRecord);
FASTArecord * fastaAllocateRecord();
int  fastaPrintRecord(FILE *ofp, FASTArecord *fRecord);
void fastaClearRecord(FASTArecord *fRecord);
void fastaDeallocateRecord(FASTArecord *fRecord);

#endif /* __FASTA_RECORD_TOOLS_HEADER__ */
