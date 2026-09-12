#pragma once
#include <stdbool.h>

// Internal, thread-local collection for a disposable, no-output geometry trial.
// No environment variable enables collection. The owner MUST discard the trial
// and perform ordinary hard-guarded construction before exposing any geometry.
// Non-wall failures retain their normal behavior.
void tok_wall_trial_begin(int movable_radial_edge);
void tok_wall_trial_capture_requested(int block, double rho);
long tok_wall_trial_end(void);
bool tok_wall_trial_record(bool fixed_radial_boundary);
bool tok_wall_trial_has_fixed_violation(void);
