/**
 * @file bssnChkptSchema.h
 * @brief Single source of truth for checkpoint (.cp) metadata and filenames.
 *
 * The field list below is visited by both a writer and a reader, so the two
 * directions cannot disagree. Add a new checkpointed field in chkpt_fields()
 * and both the CPU and GPU contexts pick it up. BH history is NOT here -- it
 * has three legacy on-disk formats and keeps its own encode/restore pair.
 */

#pragma once

#include <mpi.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

#include "TreeNode.h"
#include "dendro.h"
#include "json.hpp"

namespace bssn {

using json = nlohmann::json;

/**@brief scalar metadata carried in a .cp file. */
struct ChkptMeta {
    double tb                 = 0.0;
    double te                 = 0.0;
    double t                  = 0.0;
    double th                 = 0.0;
    unsigned int step         = 0;
    unsigned int elementOrder = 0;
    double waveletTol         = 0.0;
    double loadImbTol         = 0.0;
    unsigned int numVars      = 0;
    unsigned int activeCommSz = 0;
    double bh1[3]             = {0.0, 0.0, 0.0};
    double bh2[3]             = {0.0, 0.0, 0.0};
    double bhv1[3]            = {0.0, 0.0, 0.0};
    double bhv2[3]            = {0.0, 0.0, 0.0};
    bool bhVelValid           = false;

    bool bhMerged             = false;
    double mergeTime          = 0.0;
    unsigned int mergeStep    = 0;
    bool mergedChkptWritten   = false;

    // min det(gt) before enforcement, not yet written to _Constraints.dat
    double minDetGt           = std::numeric_limits<double>::max();

    // octree depth the .oct coordinates were written at
    unsigned int octreeDepth  = 0;

    // set by the reader: was the optional group present on disk?
    bool hasBhMerge           = false;
    bool hasMergedLatch       = false;
};

/**@brief THE field list. Add new checkpoint fields here and nowhere else.
 * req() must exist on read, opt() may be absent (older checkpoints), wonly()
 * is written but never read back. */
template <typename V>
void chkpt_fields(V& v, ChkptMeta& m) {
    v.req("DENDRO_TS_TIME_BEGIN", m.tb);
    v.req("DENDRO_TS_TIME_END", m.te);
    v.req("DENDRO_TS_ELEMENT_ORDER", m.elementOrder);

    v.req("DENDRO_TS_TIME_CURRENT", m.t);
    v.req("DENDRO_TS_STEP_CURRENT", m.step);
    v.req("DENDRO_TS_TIME_STEP_SIZE", m.th);
    v.wonly("DENDRO_TS_LAST_IO_TIME", m.t);

    v.req("DENDRO_TS_WAVELET_TOLERANCE", m.waveletTol);
    v.req("DENDRO_TS_LOAD_IMB_TOLERANCE", m.loadImbTol);
    v.req("DENDRO_TS_NUM_VARS", m.numVars);
    v.req("DENDRO_TS_ACTIVE_COMM_SZ", m.activeCommSz);

    v.req("DENDRO_BH1_X", m.bh1[0]);
    v.req("DENDRO_BH1_Y", m.bh1[1]);
    v.req("DENDRO_BH1_Z", m.bh1[2]);
    v.req("DENDRO_BH2_X", m.bh2[0]);
    v.req("DENDRO_BH2_Y", m.bh2[1]);
    v.req("DENDRO_BH2_Z", m.bh2[2]);

    v.opt("DENDRO_BH_VEL_VALID", m.bhVelValid);
    v.opt("DENDRO_BH1_VX", m.bhv1[0]);
    v.opt("DENDRO_BH1_VY", m.bhv1[1]);
    v.opt("DENDRO_BH1_VZ", m.bhv1[2]);
    v.opt("DENDRO_BH2_VX", m.bhv2[0]);
    v.opt("DENDRO_BH2_VY", m.bhv2[1]);
    v.opt("DENDRO_BH2_VZ", m.bhv2[2]);

    m.hasBhMerge = v.opt("DENDRO_BSSN_BH_MERGE", m.bhMerged);
    v.opt("DENDRO_BSSN_BH_MERGE_TIME", m.mergeTime);
    v.opt("DENDRO_BSSN_BH_MERGE_STEP", m.mergeStep);

    m.hasMergedLatch =
        v.opt("DENDRO_BSSN_MERGED_CHKPT_WRITTEN", m.mergedChkptWritten);

    v.opt("DENDRO_BSSN_MIN_DET_GT", m.minDetGt);

    v.opt("DENDRO_BSSN_OCTREE_DEPTH", m.octreeDepth);
}

/**@brief visitor that serializes into a json object. */
class ChkptWriter {
   public:
    explicit ChkptWriter(json& j) : m_j(j) {}

    template <typename T>
    bool req(const char* key, const T& val) {
        m_j[key] = val;
        return true;
    }
    template <typename T>
    bool opt(const char* key, const T& val) {
        m_j[key] = val;
        return true;
    }
    template <typename T>
    bool wonly(const char* key, const T& val) {
        m_j[key] = val;
        return true;
    }

   private:
    json& m_j;
};

/**@brief visitor that deserializes from a json object. Missing optional keys
 * leave the field at whatever the caller seeded it with. */
class ChkptReader {
   public:
    explicit ChkptReader(const json& j) : m_j(j) {}

    template <typename T>
    bool req(const char* key, T& val) {
        val = m_j.at(key).template get<T>();
        return true;
    }
    template <typename T>
    bool opt(const char* key, T& val) {
        const auto it = m_j.find(key);
        if (it == m_j.end()) return false;
        val = it->template get<T>();
        return true;
    }
    template <typename T>
    bool wonly(const char*, const T&) {
        return false;
    }

   private:
    const json& m_j;
};

inline void chkpt_write_meta(json& j, ChkptMeta& m) {
    ChkptWriter w(j);
    chkpt_fields(w, m);
}

inline void chkpt_read_meta(const json& j, ChkptMeta& m) {
    ChkptReader r(j);
    chkpt_fields(r, m);
}

// ---------------------------------------------------------------------------
// Canonical checkpoint filenames. The GPU context used to spell .cp and .oct
// differently, which made CPU and GPU checkpoints mutually unreadable; both now
// write these, and restore falls back to the legacy GPU spellings.
// ---------------------------------------------------------------------------

inline void chkpt_fname_step(char* out, size_t n, const std::string& prefix,
                             unsigned int slot) {
    snprintf(out, n, "%s_%u_step.cp", prefix.c_str(), slot);
}

inline void chkpt_fname_oct(char* out, size_t n, const std::string& prefix,
                            unsigned int slot, unsigned int rank) {
    snprintf(out, n, "%s_%u_octree_%u.oct", prefix.c_str(), slot, rank);
}

inline void chkpt_fname_var(char* out, size_t n, const std::string& prefix,
                            unsigned int slot, unsigned int rank) {
    snprintf(out, n, "%s_%u_%u.var", prefix.c_str(), slot, rank);
}

/**@brief pre-unification GPU spellings, still readable on restore. */
inline void chkpt_fname_step_legacy_gpu(char* out, size_t n,
                                        const std::string& prefix,
                                        unsigned int slot) {
    snprintf(out, n, "%s_step_%u.cp", prefix.c_str(), slot);
}

inline void chkpt_fname_oct_legacy_gpu(char* out, size_t n,
                                       const std::string& prefix,
                                       unsigned int slot, unsigned int rank) {
    snprintf(out, n, "%s_octree_%u_%u.oct", prefix.c_str(), slot, rank);
}

/**@brief sentinel naming the newest fully-written normal slot. */
inline void chkpt_fname_latest(char* out, size_t n, const std::string& prefix) {
    snprintf(out, n, "%s.latest", prefix.c_str());
}

/**@brief resolve a slot's .cp, canonical spelling first then the legacy GPU
 * one. On failure `out` keeps the canonical name, for error messages. */
inline bool chkpt_resolve_step(char* out, size_t n, const std::string& prefix,
                               unsigned int slot) {
    chkpt_fname_step(out, n, prefix, slot);
    if (std::filesystem::exists(out)) return true;

    char alt[512];
    chkpt_fname_step_legacy_gpu(alt, sizeof(alt), prefix, slot);
    if (std::filesystem::exists(alt)) {
        snprintf(out, n, "%s", alt);
        return true;
    }

    chkpt_fname_step(out, n, prefix, slot);
    return false;
}

/**@brief same, for the per-rank .oct file. */
inline bool chkpt_resolve_oct(char* out, size_t n, const std::string& prefix,
                              unsigned int slot, unsigned int rank) {
    chkpt_fname_oct(out, n, prefix, slot, rank);
    if (std::filesystem::exists(out)) return true;

    char alt[512];
    chkpt_fname_oct_legacy_gpu(alt, sizeof(alt), prefix, slot, rank);
    if (std::filesystem::exists(alt)) {
        snprintf(out, n, "%s", alt);
        return true;
    }

    chkpt_fname_oct(out, n, prefix, slot, rank);
    return false;
}

/**@brief move a restored octree from the depth it was written at to
 * m_uiMaxDepth: a deeper octree shifts every coordinate left, a shallower one
 * shifts right, levels unchanged. Fails (rank 0 says why) if an octant is finer
 * than a mesh at m_uiMaxDepth can hold. Collective over comm. */
inline bool chkpt_rescale_octree(std::vector<ot::TreeNode>& octree,
                                 unsigned int fileDepth, MPI_Comm comm) {
    if (fileDepth == m_uiMaxDepth) return true;

    int rank;
    MPI_Comm_rank(comm, &rank);

    unsigned int lmax = 0;
    for (const auto& oct : octree) lmax = std::max(lmax, oct.getLevel());
    unsigned int lmaxGlobal = 0;
    MPI_Allreduce(&lmax, &lmaxGlobal, 1, MPI_UNSIGNED, MPI_MAX, comm);

    // finest octant level remeshing allows at m_uiMaxDepth
    if (lmaxGlobal + MAXDEAPTH_LEVEL_DIFF + 1 > m_uiMaxDepth) {
        if (!rank)
            std::cout << "[BSSNCtx] : Checkpoint octree (depth " << fileDepth
                      << ") refines to level "
                      << lmaxGlobal + MAXDEAPTH_LEVEL_DIFF + 1
                      << ", finer than BSSN_MAXDEPTH " << m_uiMaxDepth
                      << " allows" << std::endl;
        return false;
    }

    const bool deeper = fileDepth < m_uiMaxDepth;
    const unsigned int k =
        deeper ? m_uiMaxDepth - fileDepth : fileDepth - m_uiMaxDepth;
    for (auto& oct : octree) {
        const unsigned int x = deeper ? oct.getX() << k : oct.getX() >> k;
        const unsigned int y = deeper ? oct.getY() << k : oct.getY() >> k;
        const unsigned int z = deeper ? oct.getZ() << k : oct.getZ() >> k;
        ot::TreeNode moved(x, y, z, 0, m_uiDim, m_uiMaxDepth);
        moved.setFlag(oct.getFlag());
        moved.setWeight(oct.getWeight());
        oct = moved;
    }

    if (!rank)
        std::cout << "[BSSNCtx] : Checkpoint octree moved from depth "
                  << fileDepth << " to BSSN_MAXDEPTH " << m_uiMaxDepth
                  << std::endl;
    return true;
}

/**@brief bumped only when the on-disk layout changes incompatibly. */
constexpr int BSSN_CHKPT_FORMAT_VERSION = 1;

/**@brief publish the sentinel via .tmp + atomic rename, so it only ever names
 * a checkpoint that finished writing. Rank 0 only; no-op on other ranks. */
inline void chkpt_publish_latest(const std::string& prefix, unsigned int slot,
                                 unsigned int step, unsigned int rank) {
    if (rank) return;

    char finalf[512];
    char tmpf[512];
    chkpt_fname_latest(finalf, sizeof(finalf), prefix);
    snprintf(tmpf, sizeof(tmpf), "%s.latest.tmp", prefix.c_str());

    json latest;
    latest["format_version"] = BSSN_CHKPT_FORMAT_VERSION;
    latest["slot"]           = slot;
    latest["step"]           = step;

    {
        std::ofstream lf(tmpf);
        if (!lf) return;  // sentinel is an optimization; restore falls back
        lf << latest << std::endl;
    }

    std::error_code ec;
    std::filesystem::rename(tmpf, finalf, ec);  // atomic publish
}

/**@brief read the sentinel's slot. Returns false when absent/unreadable/from a
 * newer format, which sends restore back to the legacy 0/1 scan. */
inline bool chkpt_read_latest(const std::string& prefix, unsigned int& slot) {
    char fname[512];
    chkpt_fname_latest(fname, sizeof(fname), prefix);
    if (!std::filesystem::exists(fname)) return false;

    std::ifstream lf(fname);
    if (!lf) return false;

    json latest;
    try {
        lf >> latest;
    } catch (...) {
        return false;
    }

    if (latest.value("format_version", -1) > BSSN_CHKPT_FORMAT_VERSION)
        return false;

    const int s = latest.value("slot", -1);
    if (s < 0) return false;

    slot = (unsigned int)s;
    return true;
}

/**@brief true if this rank's .oct and .var for the slot are as long as their
 * own headers say, i.e. not the short files a failed write leaves behind. */
inline bool chkpt_rank_files_complete(const std::string& prefix,
                                      unsigned int slot, unsigned int rank,
                                      unsigned int numVars) {
    char f[512];
    std::error_code ec;

    unsigned int num = 0;
    chkpt_resolve_oct(f, sizeof(f), prefix, slot, rank);
    FILE* fp = fopen(f, "rb");
    if (!fp) return false;
    const bool octHdr = fread(&num, sizeof(num), 1, fp) == 1;
    fclose(fp);
    if (!octHdr) return false;
    const uintmax_t octSz = std::filesystem::file_size(f, ec);
    if (ec ||
        octSz != sizeof(unsigned int) + (uintmax_t)num * sizeof(ot::TreeNode))
        return false;

    // header: total nodes, local begin, local end
    unsigned int hdr[3];
    chkpt_fname_var(f, sizeof(f), prefix, slot, rank);
    fp = fopen(f, "rb");
    if (!fp) return false;
    const bool varHdr = fread(hdr, sizeof(unsigned int), 3, fp) == 3;
    fclose(fp);
    if (!varHdr || hdr[2] < hdr[1]) return false;
    const uintmax_t varSz = std::filesystem::file_size(f, ec);
    const uintmax_t body =
        hdr[0] > 0 ? (uintmax_t)numVars * (hdr[2] - hdr[1]) * sizeof(double)
                   : 0;
    return !ec && varSz == 3 * sizeof(unsigned int) + body;
}

/**@brief true if every rank that wrote the slot left complete files. Collective
 * over comm; the writer count comes from the slot's .cp. */
inline bool chkpt_slot_complete(const std::string& prefix, unsigned int slot,
                                unsigned int numVars, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    unsigned int writers = 0;
    if (!rank) {
        char f[512];
        if (chkpt_resolve_step(f, sizeof(f), prefix, slot)) {
            std::ifstream in(f);
            json cp;
            try {
                in >> cp;
                writers = cp.value("DENDRO_TS_ACTIVE_COMM_SZ", 0u);
            } catch (...) {
                writers = 0;
            }
        }
    }
    MPI_Bcast(&writers, 1, MPI_UNSIGNED, 0, comm);
    if (writers == 0) return false;

    int ok    = ((unsigned int)rank >= writers) ||
                chkpt_rank_files_complete(prefix, slot, rank, numVars);
    int allOk = 0;
    MPI_Allreduce(&ok, &allOk, 1, MPI_INT, MPI_MIN, comm);
    return allOk == 1;
}

}  // namespace bssn
