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

/** Record a wall violation, saying whether the offending node was itself
 * outside the vessel (as opposed to an edge between two inside nodes bulging
 * out). Only the first kind, on a non-movable radial boundary, means the
 * plasma boundary leaves the machine. */
bool tok_wall_trial_record_scope(bool fixed_radial_boundary, bool node_outside);

/** True when a node of the NON-movable radial boundary was itself outside the
 * vessel: the separatrix row leaves the machine, so the configuration is
 * inner-wall-limited rather than diverted. No boundary adjustment can change
 * that -- the offending row is a declared join, never a movable edge. */
bool tok_wall_trial_has_fixed_node_outside(void);
