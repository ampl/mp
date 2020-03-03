#include <omp.h>

 static omp_lock_t Locks[5];

 void
INIT_DTOA_LOCKS(void)
{
	int i;
	static int need_init = 1;

	if (need_init) {
		for(i = 0; i < 5; ++i)
			omp_init_lock(&Locks[i]);
		need_init = 0;
		}
	}

 void
ACQUIRE_DTOA_LOCK(unsigned int n)
{
	if (n < 5)
		omp_set_lock(&Locks[n]);
	}

 void
FREE_DTOA_LOCK(unsigned int n)
{
	if (n < 5)
		omp_unset_lock(&Locks[n]);
	}
