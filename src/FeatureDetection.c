#include "FeatureDetection.h"
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

#if defined(_WIN32)
    #include <windows.h>
    #include <intrin.h>
#elif defined(__linux__) || defined(__unix__) || defined(__posix__)
    #include <unistd.h>
    #include <dirent.h>
    #include <cpuid.h>
#endif

struct system_features_t features = {0};

#if defined(__linux__) || defined(__unix__) || defined(__posix__)

static int count_directories(const char *path) {
    DIR *dir;
    struct dirent *entry;
    int count = 0;

    dir = opendir(path);
    if (!dir) {
        perror("Error opening directory.");
        return -1;
    }

    while ((entry = readdir(dir)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
            continue;
        if (entry->d_type == DT_DIR)
            ++count;
    }

    closedir(dir);
    return count;
}

static int query_file_size(const char *path) {
    FILE *fp;
    char buff[64];
    int size = 0;

    fp = fopen(path, "r");
    if (!fp) return -1;

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
    char buff[64];
    int shared = -1;

    fp = fopen(path, "r");
    if (!fp) return -1;

    if (fgets(buff, sizeof(buff), fp)) {
        shared = strchr(buff, '-') ? 1 : 0;
    }

    fclose(fp);
    return shared;
}

static int query_file_type(const char *path) {
    FILE *fp;
    char buff[64];

    fp = fopen(path, "r");
    if (!fp) return -1;

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
    char buff[64];
    int level = -1;

    fp = fopen(path, "r");
    if (!fp) return -1;

    if (fgets(buff, sizeof(buff), fp)) {
        sscanf(buff, "%d", &level);
    }

    fclose(fp);
    return level;
}

#endif

static void query_logical_cores(struct system_features_t *f) {
#if defined(_WIN32)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    f->logical_cores = (uint8_t)si.dwNumberOfProcessors;
#elif defined(__linux__) || defined(__unix__) || defined(__posix__)
    f->logical_cores = sysconf(_SC_NPROCESSORS_ONLN);
#else
    f->logical_cores = 0;
#endif
}

static void query_physical_cores(struct system_features_t *f) {
#if defined(_WIN32)
    DWORD len = 0;
    GetLogicalProcessorInformation(NULL, &len);
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION info = malloc(len);
    if (!info || !GetLogicalProcessorInformation(info, &len)) {
        f->physical_cores = f->logical_cores;
    } else {
        DWORD count = len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
        DWORD physical = 0;
        for (DWORD i = 0; i < count; ++i) {
            if (info[i].Relationship == RelationProcessorCore)
                ++physical;
        }
        f->physical_cores = (uint8_t)physical;
    }
    free(info);
#elif defined(__linux__)
    FILE *fp = fopen("/proc/cpuinfo", "r");
    if (!fp) {
        f->physical_cores = sysconf(_SC_NPROCESSORS_ONLN);
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
    f->physical_cores = (cores > 0) ? cores : sysconf(_SC_NPROCESSORS_ONLN);
#else
    f->physical_cores = 0;
#endif
}

static void query_cache_sizes(struct system_features_t *f) {
    memset(f->l, 0, sizeof(f->l));
#if defined(_WIN32)
    DWORD len = 0;
    GetLogicalProcessorInformation(NULL, &len);
    PSYSTEM_LOGICAL_PROCESSOR_INFORMATION info = malloc(len);
    if (!info || !GetLogicalProcessorInformation(info, &len)) {
        f->line_size = 0;
        free(info);
        return;
    }

    DWORD count = len / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
    for (DWORD i = 0; i < count; ++i) {
        if (info[i].Relationship == RelationCache) {
            CACHE_DESCRIPTOR *cd = &info[i].Cache;
            int level = cd->Level;
            f->l[level] = cd->Size;
            f->line_size = cd->LineSize;
            if (level == 2)
                f->l2_shared = (cd->Associativity > 1);
        }
    }
    free(info);
#elif defined(__linux__)
    char *base = "/sys/devices/system/cpu/cpu0/cache/";
    int dc = count_directories(base);
    for (int i = 0; i < dc; ++i) {
        char path[256];
        size_t ps = sizeof(path);
        snprintf(path, ps, "%sindex%d/size", base, i);
        int size = query_file_size(path);
        snprintf(path, ps, "%sindex%d/shared_cpu_list", base, i);
        int shared = query_file_shared(path);
        snprintf(path, ps, "%sindex%d/type", base, i);
        int data = query_file_type(path);
        snprintf(path, ps, "%sindex%d/level", base, i);
        int level = query_file_gen(path);
        snprintf(path, ps, "%sindex%d/coherency_line_size", base, i);
        int line = query_file_gen(path);

        if (size == -1 || shared == -1 || data == -1 || level == -1)
            continue;

        if (data)
            f->l[level] = size;

        if (shared && level == 2)
            f->l2_shared = 1;
        else if (level == 2)
            f->l2_shared = 0;

        f->line_size = line;
    }
#endif
}

static void query_cpu_features(struct system_features_t *f) {
#if defined(_WIN32) || defined(__x86_64__) || defined(_M_X64)
    int regs[4];
#if defined(_WIN32)
    __cpuid(regs, 1);
#else
    __cpuid(1, regs[0], regs[1], regs[2], regs[3]);
#endif
    uint32_t ecx = regs[2];

    f->sse = (ecx & (1u << 19)) != 0;
    int has_avx = (ecx & (1u << 28)) != 0;   
    int has_fma = (ecx & (1u << 12)) != 0;

#if defined(_WIN32)
    uint64_t xcr0 = _xgetbv(0);
#else
    uint32_t xcr0_lo, xcr0_hi;
    __asm__ volatile(".byte 0x0f, 0x01, 0xd0"
                     : "=a"(xcr0_lo), "=d"(xcr0_hi)
                     : "c"(0));
    uint64_t xcr0 = ((uint64_t)xcr0_hi << 32) | xcr0_lo;
#endif

    int os_avx = ((xcr0 & 0x6) == 0x6);
    f->avx = has_avx && os_avx;
    f->fma3 = f->avx && has_fma;
#endif
}

void query_features(struct system_features_t *f) {
    memset(f, 0, sizeof(*f));
    query_logical_cores(f);
    query_physical_cores(f);
    query_cache_sizes(f);
    query_cpu_features(f);
    f->initialized = true;
    if (f->avx) {
      f->instructions = AVX;
    } else if (f->sse) {
      f->instructions = SSE;
    } else {
      f->instructions = FB;
    }
}
