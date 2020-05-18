#ifndef COMMONS_INCLUDED
#define COMMONS_INCLUDED
// Functions that may be used by different .cpp files

void matrix_to_tensor(const arma::imat& Z, arma::icube& Z_tensor);
void vb_matrix_to_tensor(const arma::imat& Z, arma::cube& Z_tensor);


#endif