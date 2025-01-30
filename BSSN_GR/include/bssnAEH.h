#pragma once
#include "aeh.h"
#include "bssnCtx.h"

namespace bssnaeh {

extern std::vector<DendroScalar> hh[4];

void initialize_aeh();

void perform_aeh_step(bssn::BSSNCtx* const bssnCtx, const int rank);

}  // namespace bssnaeh
