#include <stdlib.h>

void stime(char * buffer, int size);

#ifdef _WIN32
int mkstemp(char *tmpl);
#endif

void mkdir_p(const char *path);

