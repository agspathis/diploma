#include <fstream>

#include "utilities.h"

float rand_fb (float b)
{
    return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/b));
}

int print_aabb(aabb aabb)
{
    printf("min:%f %f %f\n", aabb.min.getX(), aabb.min.getY(), aabb.min.getZ());
    printf("max:%f %f %f\n", aabb.max.getX(), aabb.max.getY(), aabb.max.getZ());
    return 0;
}

int print_long_array(long* array, long count)
{
    for (long index=0; index<count; index++) {
	printf("%d ", array[index]);
	if (index % 16) continue;
	printf("\n");
    }
    printf("\n");
    return 0;
}
