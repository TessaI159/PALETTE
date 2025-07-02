#include "FeatureDetection.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__) || defined(__unix__) || defined(__posix__)
#include <dirent.h>
#include <unistd.h>
#endif

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <cpuid.h>
#endif
#endif

/* TODO (Tess): In the Linux cache detection, check for errors for each of the
   mini query functions*/

#if defined(__linux__) || defined(__unix__) || defined(__poxis__)

static int count_directories(const char *path) {
	DIR	      *dir;
	struct dirent *entry;
	int	       count = 0;

	dir = opendir(path);
	if (!dir) {
		perror("Error opening directory.");
		return -1;
	}

	while ((entry = readdir(dir)) != NULL) {
		if (strcmp(entry->d_name, ".") == 0 ||
		    strcmp(entry->d_name, "..") == 0) {
			continue;
		}
		if (entry->d_type == DT_DIR) {
			++count;
		}
	}
	closedir(dir);
	return count;
}

static int query_file_size(const char *path) {
	FILE *fp;
	char  buff[64];
	fp = fopen(path, "r");
	if (!fp) {
		printf("Could not open %s\n", path);
		return -1;
	}

	int size = 0;
	if (fgets(buff, sizeof(buff), fp)) {
		if (strchr(buff, 'K')) {
			sscanf(buff, "%dK", &size);
			size *= 1024;
		} else if (strchr(buff, 'M')) {
			sscanf(buff, "%dM", &size);
			size *= 1024 * 1024;
		}
	}
	fclose(fp);
	return size;
}

static int query_file_shared(const char *path) {
	FILE *fp;
	char  buff[64];
	fp = fopen(path, "r");
	if (!fp) {
		printf("Could not open %s\n", path);
		return -1;
	}
	int shared = -1;
	if (fgets(buff, sizeof(buff), fp)) {
		if (strchr(buff, '-')) {
			shared = 1;
		} else {
			shared = 0;
		}
	}
	fclose(fp);
	return shared;
}

static int query_file_type(const char *path) {
	FILE *fp;
	char  buff[64];
	fp = fopen(path, "r");
	if (!fp) {
		printf("Could not open %s\n", path);
		return -1;
	}

	if (fgets(buff, sizeof(buff), fp)) {
		if (strstr(buff, "Data") || strstr(buff, "Unified")) {
			fclose(fp);
			return 1;
		}
	}
	fclose(fp);
	return -1;
}

static int query_file_gen(const char *path) {
	FILE *fp;
	char  buff[64];
	fp = fopen(path, "r");
	int level;
	if (!fp) {
		printf("Could not open %s\n", path);
		return -1;
	}

	if (fgets(buff, sizeof(buff), fp)) {
		sscanf(buff, "%d", &level);
	}
	fclose(fp);
	return level;
}

#endif

static void query_cache_sizes(struct system_features_t *features) {
	features->l[0] = features->l[1] = features->l[2] = 0;
#if defined(__linux__)
	char *base = "/sys/devices/system/cpu/cpu0/cache/";
	int   dc   = count_directories(base);

	for (int i = 0; i < dc; ++i) {
		char path[256];
		snprintf(path, sizeof(path), "%sindex%d/size", base, i);
		int size = query_file_size(path);
		snprintf(path, sizeof(path), "%sindex%d/shared_cpu_list", base,
			 i);
		int shared = query_file_shared(path);
		snprintf(path, sizeof(path), "%sindex%d/type", base, i);
		int data = query_file_type(path);
		snprintf(path, sizeof(path), "%sindex%d/level", base, i);
		int level = query_file_gen(path);
		snprintf(path, sizeof(path), "%sindex%d/coherency_line_size",
			 base, i);
		int line = query_file_gen(path);

		if (size == -1 || shared == -1 || data == -1 || level == -1)
			continue;

		if (data)
			features->l[level - 1] = size;

		if (shared && level == 2)
			features->l2_shared = 1;
		else
			features->l2_shared = 0;
		features->line_size = line;
	}
#elif defined(__x86_64__) || defined(_M_X64)
	for (int i = 0; i < 10; ++i) {
		uint32_t eax, ebx, ecx, edx;
		__cpuid_count(4, i, eax, ebx, ecx, edx);

		uint32_t cache_type = eax & 0x1F;
		if (cache_type == 0)
			break;

		uint32_t level	    = (eax >> 5) & 0x7;
		uint32_t ways	    = ((ebx >> 22) & 0x3FF) + 1;
		uint32_t partitions = ((ebx >> 12) & 0x3FF) + 1;
		uint32_t line_size  = (ebx & 0xFFF) + 1;
		uint32_t sets	    = ecx + 1;
		uint32_t size_bytes = ways * partitions * line_size * sets;

		if (level == 1)
			features->l[0] += size_bytes;
		else if (level == 2)
			features->l[1] += size_bytes;
		else if (level == 3)
			features->l[2] += size_bytes;
	}

#endif
}

static void query_logical_cores(struct system_features_t *features) {
#if defined(_WIN32)
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	features->logical_cores = sysinfo.dwNumberOfProcessors;
#elif defined(__linux__) || defined(__unix__) || defined(__posix__)
	features->logical_cores = sysconf(_SC_NPROCESSORS_ONLN);
#else
	features->logical_cores = -1;
#endif
}

static void query_physical_cores(struct system_features_t *features) {
#if defined(__linux__)
	FILE *fp = fopen("/proc/cpuinfo", "r");
	if (!fp) {
		features->physical_cores = sysconf(_SC_NPROCESSORS_ONLN);
		return;
	}
	int  physical_ids[256], core_ids[256], cores = 0;
	int  physical_id = -1, core_id = -1;
	char line[256];
	while (fgets(line, sizeof(line), fp)) {
		if (strncmp(line, "physical id", 11) == 0)
			sscanf(line, "physical id : %d", &physical_id);
		else if (strncmp(line, "core id", 8) == 0) {
			sscanf(line, "core id : %d", &core_id);
			if (physical_id >= 0 && core_id >= 0) {
				int unique = 1;
				for (int i = 0; i < cores; ++i) {
					if (physical_ids[i] == physical_id &&
					    core_ids[i] == core_id) {
						unique = 0;
						break;
					}
				}
				if (unique && cores < 256) {
					physical_ids[cores] = physical_id;
					core_ids[cores]	    = core_id;
					++cores;
				}
				physical_id = core_id = -1;
			}
		}
	}
	fclose(fp);
	features->physical_cores =
	    (cores > 0) ? cores : sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_WIN32)
	DWORD len = 0;
	GetLogicalProcessorInformationEx(RelationProcessorCore, NULL, &len);
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *buffer =
	    (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *)malloc(len);
	if (!buffer) {
		features->physical_cores = -1;
		return;
	}
	if (!GetLogicalProcessorInformationEx(RelationProcessorCore, buffer,
					      &len)) {
		free(buffer);
		features->physical_cores = -1;
		return;
	}
	int   count = 0;
	char *ptr   = (char *)buffer;
	while (ptr < (char *)buffer + len) {
		SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *item =
		    (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *)ptr;
		if (item->Relationship == RelationProcessorCore)
			++count;
		ptr += item->Size;
	}
	free(buffer);
	features->physical_cores = count;
#else
	features->physical_cores = -1;
#endif
}

static void query_cpu_features(struct system_features_t *features) {
	features->avx = features->avx2 = features->fma3 = 0;
#if defined(__x86_64__) || defined(_M_X64)
	uint32_t eax, ebx, ecx, edx;
	eax = ebx = ecx = edx = 0;
#if defined(_MSC_VER)
	int cpu_info[4];
	__cpuid(cpu_info, 1);
	ecx = cpu_info[2];
#else
	__cpuid(1, eax, ebx, ecx, edx);
#endif
	int has_avx = (ecx & (1u << 28)) != 0;
	int has_fma = (ecx & (1u << 12)) != 0;

	uint32_t xcr0_lo = 0, xcr0_hi = 0;
#if defined(_MSC_VER)
	xcr0_lo = (uint32_t)_xgetbv(0);
#else
	__asm__ volatile(".byte 0x0f, 0x01, 0xd0"
			 : "=a"(xcr0_lo), "=d"(xcr0_hi)
			 : "c"(0));
#endif
	int os_supports_avx = (xcr0_lo & 0x6) == 0x6;

	features->avx  = has_avx && os_supports_avx;
	features->fma3 = features->avx && has_fma;

	if (features->avx) {
#if defined(_MSC_VER)
		__cpuid(cpu_info, 7);
		ebx = cpu_info[1];
#else
		__cpuid_count(7, 0, eax, ebx, ecx, edx);
#endif
		features->avx2 = (ebx & (1 << 5)) != 0;
	}
#endif
}

void query_features(struct system_features_t *features) {
	memset(features, 0, sizeof(struct system_features_t));
	query_logical_cores(features);
	query_physical_cores(features);
	query_cache_sizes(features);
	query_cpu_features(features);
}
