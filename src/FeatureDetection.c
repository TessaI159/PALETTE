#include "FeatureDetection.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
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

void query_cache_sizes(struct system_features_t* features) {
	features->l1 = features->l2 = features->l3 = 0;
#if defined(__x86_64__) || defined(_M_X64)
	for (int i = 0; i < 10; ++i) {
		uint32_t eax, ebx, ecx, edx;
		__cpuid_count(4, i, eax, ebx, ecx, edx);

		uint32_t cache_type = eax & 0x1F;
		if (cache_type == 0)
			break;

		uint32_t level = (eax >> 5) & 0x7;
		uint32_t ways = ((ebx >> 22) & 0x3FF) + 1;
		uint32_t partitions = ((ebx >> 12) & 0x3FF) + 1;
		uint32_t line_size = (ebx & 0xFFF) + 1;
		uint32_t sets = ecx + 1;
		uint32_t size_bytes = ways * partitions * line_size * sets;

		if (level == 1)
			features->l1 += size_bytes;
		else if (level == 2)
			features->l2 += size_bytes;
		else if (level == 3)
			features->l3 += size_bytes;
	}
#elif defined(__linux__)
	FILE* fp;
	char path[256], buff[32];
	for (int i = 0; i < 3; ++i) {
		snprintf(path, sizeof(path),
			 "/sys/devices/system/cpu/cpu0/cache/index%d/size", i);
		fp = fopen(path, "r");
		if (!fp)
			continue;

		if (fgets(buff, sizeof(buff), fp)) {
			int size = 0;
			if (strchr(buff, 'K')) {
				sscanf(buff, "%dK", &size);
				size *= 1024;
			} else if (strchr(buff, 'M')) {
				sscanf(buff, "%dM", &size);
				size *= 1024 * 1024;
			}
			if (i == 0)
				features->l1 += size;
			else if (i == 1)
				features->l2 += size;
			else if (i == 2)
				features->l3 += size;
		}
		fclose(fp);
	}

#endif
}

void query_logical_cores(struct system_features_t* features) {
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

void query_physical_cores(struct system_features_t* features) {
#if defined(__linux__)
	FILE* fp = fopen("/proc/cpuinfo", "r");
	if (!fp) {
		features->physical_cores = sysconf(_SC_NPROCESSORS_ONLN);
		return;
	}
	int physical_ids[256], core_ids[256], cores = 0;
	int physical_id = -1, core_id = -1;
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
					core_ids[cores] = core_id;
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
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* buffer =
	    (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)malloc(len);
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
	int count = 0;
	char* ptr = (char*)buffer;
	while (ptr < (char*)buffer + len) {
		SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* item =
		    (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)ptr;
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

void query_cpu_features(struct system_features_t* features) {
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

	features->avx = has_avx && os_supports_avx;
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
