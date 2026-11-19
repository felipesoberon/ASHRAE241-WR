#ifndef ECAI_CPP2_PAIRED_OUTPUT_H
#define ECAI_CPP2_PAIRED_OUTPUT_H

#include "cli.h"

namespace ecai_cpp2::paired_output {
// Writes to a new directory. A manifest is published only after all data closes.
void run_and_write(const cli::ResolvedRun& run, const cli::CliOptions& options,
                   const char* protected_legacy_tree);
}
#endif
