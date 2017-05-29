
#define MAX_FILENAME 1024

enum bamdb_convert_to {
	BAMDB_CONVERT_TO_TEXT,
	BAMDB_CONVERT_TO_SQLITE,
	BAMDB_CONVERT_TO_LMDB
};

typedef struct bamdb_args {
	enum bamdb_convert_to convert_to;
	char input_file_name[MAX_FILENAME];
	char *index_file_name;
        char *output_file_name;
	char *bx;
} bam_args_t;

int generate_index_file(char *input_file_name);
