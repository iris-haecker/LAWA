namespace lawa {

template <typename T, typename Basis2D, typename Preconditioner>
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::AdaptiveHelmholtzOperator2D
                                                         (const Basis2D &_basis, T _c,
                                                          const Preconditioner &_Prec,
                                                          T thresh, int NumOfCols, int NumOfRows)
: basis(_basis), c(_c), Prec(_Prec),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  op_identity_x(basis.first), op_identity_y(basis.second),
  op_laplace_x(basis.first), op_laplace_y(basis.second),
  Prec1D(),
  data_identity_x(op_identity_x, Prec1D, compression_1d_x, thresh, NumOfRows, NumOfCols),
  data_identity_y(op_identity_y, Prec1D, compression_1d_y, thresh, NumOfRows, NumOfCols),
  data_laplace_x(op_laplace_x,   Prec1D, compression_1d_x, thresh, NumOfRows, NumOfCols),
  data_laplace_y(op_laplace_y,   Prec1D, compression_1d_y, thresh, NumOfRows, NumOfCols),
  P_data()
{

}

template <typename T, typename Basis2D, typename Preconditioner>
T
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::operator()(const Index2D &row_index,
                                                                    const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    if (!flens::IsSame<NoPreconditioner2D, Preconditioner>::value) {
        const_coeff_it it_P_end       = P_data.end();
        const_coeff_it it_row_index   = P_data.find(row_index);
        if (it_row_index != it_P_end) {
            prec *= (*it_row_index).second;
        }
        else {
            T tmp;
            if (flens::IsSame<DiagonalHelmholtzPreconditioner2D, Preconditioner>::value) {
                T prec_dd_x = data_laplace_x(row_index.index1,row_index.index1);
                T prec_id_x = data_identity_x(row_index.index1,row_index.index1);
                T prec_id_y = data_identity_y(row_index.index2,row_index.index2);
                T prec_dd_y = data_laplace_y(row_index.index2,row_index.index2);
                tmp = 1./std::sqrt(fabs(prec_dd_x*prec_id_y + prec_id_x*prec_dd_y
                                        + c*prec_id_x*prec_id_y ));
            }
            else {
                tmp = Prec(row_index);
            }
            P_data[row_index] = tmp;
            prec *= tmp;
        }
        it_P_end       = P_data.end();
        const_coeff_it it_col_index   = P_data.find(col_index);
        if (it_col_index != it_P_end) {
            prec *= (*it_col_index).second;
        }
        else {
            T tmp;
            if (flens::IsSame<DiagonalHelmholtzPreconditioner2D, Preconditioner>::value) {
                T prec_dd_x = data_laplace_x(col_index.index1,col_index.index1);
                T prec_id_x = data_identity_x(col_index.index1,col_index.index1);
                T prec_id_y = data_identity_y(col_index.index2,col_index.index2);
                T prec_dd_y = data_laplace_y(col_index.index2,col_index.index2);
                tmp = 1./std::sqrt(fabs(prec_dd_x*prec_id_y + prec_id_x*prec_dd_y
                                        + c*prec_id_x*prec_id_y ));
            }
            else {
                tmp = Prec(col_index);
            }
            P_data[col_index] = tmp;
            prec *= tmp;
        }
    }

    T dd_x = data_laplace_x(row_index.index1,col_index.index1);
    T id_x = data_identity_x(row_index.index1,col_index.index1);

    T id_y, dd_y;
    if ( (flens::IsSame<typename Basis2D::FirstBasisType, Basis<T,Primal,R,CDF> >::value) &&
         (flens::IsSame<typename Basis2D::SecondBasisType, Basis<T,Primal,R,CDF> >::value)    ) {
        if (row_index.index2.xtype==XWavelet && col_index.index2.xtype==XWavelet) {
            id_y = data_identity_x(row_index.index2,col_index.index2);
            dd_y = data_laplace_x(row_index.index2,col_index.index2);
        }
        else {
            id_y = data_identity_y(row_index.index2,col_index.index2);
            dd_y = data_laplace_y(row_index.index2,col_index.index2);
        }
    }
    else {
        dd_y = data_laplace_y(row_index.index2,col_index.index2);
        id_y = data_identity_y(row_index.index2,col_index.index2);
    }

    return prec*(dd_x*id_y + id_x*dd_y + c*id_x*id_y  );
}

template <typename T, typename Basis2D, typename Preconditioner>
T
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    if (!flens::IsSame<NoPreconditioner2D, Preconditioner>::value) {
        const_coeff_it it_P_end       = P_data.end();
        const_coeff_it it_index   = P_data.find(index);
        if (it_index != it_P_end) {
            prec *= (*it_index).second;
        }
        else {
            T tmp = Prec(index);
            P_data[index] = tmp;
            prec *= tmp;
        }
    }
    return prec;
}

template <typename T, typename Basis2D, typename Preconditioner>
void
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::clear()
{
    data_laplace_x.clear();
    data_identity_x.clear();
    data_laplace_y.clear();
    data_identity_y.clear();
}


}   // namespace lawa

