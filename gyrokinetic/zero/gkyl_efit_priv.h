#include <gkyl_efit.h>
#include <gkyl_dg_basis_ops.h>

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics);

int 
find_xpts(gkyl_efit* up, double *Rxpt, double *Zxpt);

int 
find_xpts_cubic(gkyl_efit* up, double *Rxpt, double *Zxpt);

// Strips a filepath of everything through the last forward slash
// and the extension to give just the name
void
get_stripped_filename(const char *filepath, char *out_buffer);
