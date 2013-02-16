/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <list>

namespace lawa {

template <typename T, typename Basis>
void
getSingularPoints(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff, DenseVector<Array<T> > &sing_pts)
{

    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;
    std::list<T> temp;
    for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
        
        DenseVector<Array<T> > bf_singpts = basis.generator((*it).first.xtype).singularSupport((*it).first.j, (*it).first.k);

        for (int i = bf_singpts.firstIndex(); i <= bf_singpts.lastIndex(); ++i) {
            temp.push_back(bf_singpts(i));
        }
    }
    temp.sort(); temp.unique();
    sing_pts.engine().resize((int)temp.size());
    int i = 1;
    for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it ) {
        sing_pts(i) = *it; ++i;
    }
}

template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     const Preconditioner &P, T (*u)(T), const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());

    DenseVector<Array<T> > sing_pts;
    getSingularPoints(basis, coeff, sing_pts);

    for (int i=sing_pts.firstIndex(); i<=sing_pts.lastIndex(); ++i) {
        T x = sing_pts(i);
        T appr = 0.0;
        T exact= u(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j, k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            
            appr += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);

        }
        plotfile << x << " " << exact << " " << appr  << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis, typename Preconditioner>
void
plot(const Basis &basis, const Coefficients<Lexicographical,T,Index1D> coeff,
     const Preconditioner &P, T (*u)(T), T (*du)(T), T a, T b, T h, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D >::const_iterator coeff_it;

    std::ofstream plotfile(filename);
    for (T x=a; x<=b; x+=h) {
        T appr=0., d_appr = 0.0;
        T exact= u(x);
        T d_exact= du(x);
        for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
            int j = (*it).first.j, k = (*it).first.k;
            T coeff = (*it).second, prec = P((*it).first);
            
            appr   += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,0);
            d_appr += prec * coeff * basis.generator((*it).first.xtype)(x,j,k,1);

        }
        plotfile << x << " " << exact << " " << d_exact << " " << appr << " " << d_appr << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       const Preconditioner &P, T (*u)(T,T), T a1, T b1, T a2, T b2, T h, const char* filename)
{

    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());

    for (T x=a1; x<=b1; x+=h) {
        for (T y=a2; y<=b2; y+=h) {
            T appr = 0.0;
            T exact= u(x,y);
            for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype;
                XType xtype_y = (*it).first.index2.xtype;
                int j_x = (*it).first.index1.j, k_x = (*it).first.index1.k;
                int j_y = (*it).first.index2.j, k_y = (*it).first.index2.k;

                T coeff = (*it).second, prec = P((*it).first);
                
                appr += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,0);

            }
            plotfile << x << " " << y << " " << exact << " " << appr  << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}

template <typename T, typename Basis2D, typename Preconditioner>
void
plot2D(const Basis2D &basis, const Coefficients<Lexicographical,T,Index2D> coeff,
       const Preconditioner &P, T (*u)(T,T), T (*dy_u)(T,T), T a1, T b1, T a2, T b2, 
       T h1, T h2, const char* filename)
{

    typedef typename Coefficients<Lexicographical,T,Index2D >::const_iterator coeff_it;

    std::stringstream PlotFileName;
    PlotFileName << filename << ".dat";
    std::ofstream plotfile(PlotFileName.str().c_str());

    for (T x=a1; x<=b1; x+=h1) {
        for (T y=a2; y<=b2; y+=h2) {
            T appr = 0.0;
            T dy_appr = 0.0;
            T exact= u(x,y);
            T dy_exact = dy_u(x,y);
            for (coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype;
                XType xtype_y = (*it).first.index2.xtype;
                int j_x = (*it).first.index1.j, k_x = (*it).first.index1.k;
                int j_y = (*it).first.index2.j, k_y = (*it).first.index2.k;

                T coeff = (*it).second, prec = P((*it).first);
                
                appr    += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,0);
                dy_appr += prec * coeff * basis.first.generator(xtype_x)(x,j_x,k_x,0) * basis.second.generator(xtype_y)(y,j_y,k_y,1);

            }
            plotfile << x << " " << y << " " << exact << " " << appr  << " "
                     << dy_exact << " " << dy_appr << std::endl;
        }
        plotfile << std::endl;
    }
    plotfile.close();
}

template <typename T, DomainType Domain, Construction Cons>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff,
          const Basis<T,Primal,Domain,Cons> &basis,
          int j0, int J,
          const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator const_it;
    if (coeff.size() == 0) {
        return;
    }

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    gps << "reset" << std::endl;
    //gps << "set terminal postscript eps enh color; set output '" << filename << ".eps'" << std::endl;
    gps << "set terminal png; set output '" << filename << ".png'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    double width  = basis.mra.rangeI(j0).length();
    double height = (J+2);

    for (int j=j0-1; j<=std::min(j0+J,7); ++j) {
        const double fromY = j;
        const double toY   = j+1;

        int numK = (j<j0) ? basis.mra.rangeI(j0).length()
                          : basis.rangeJ(j).length();

        for (int k=1; k<=numK; ++k) {
            const double fromX = width*double(k-1)/numK;
            const double toX = width*double(k)/numK;

            gps << "set object rectangle from "
                << fromX << ", " << fromY
                << " to "
                << toX << "," << toY << " fc rgb 'black' "
                << " linewidth " << 0.1 << " fillstyle empty " << std::endl;
        }
    }

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j = (*it).second.j;
        int k = (*it).second.k;

        const double fromY = ((*it).second.xtype == XBSpline)
                           ? j0-1
                           : j;
        const double toY   = fromY + 1;

        const int numK = ((*it).second.xtype == XBSpline)
                       ? basis.mra.rangeI(j).length()
                       : basis.rangeJ(j).length();

        const double fromX = width*double(k-1)/numK;
        const double toX = width*double(k)/numK;
        
        gps << "set object rectangle from "
            << fromX << ", " << fromY
            << " to "
            << toX << "," << toY << " fc rgb 'black' "
            << " linewidth " << 0.1 << " fillstyle solid " << std::endl;
    }

    gps << "set xrange[-0.5:" << (width+0.5) << "]" << std::endl;
    gps << "set yrange[" << j0-1.5 << ":" << (j0+J+1.5) << "]" << std::endl;

    //gps << "set xtics 1" << endl;
    gps << "set ytics ('" << j0 << "' " << j0-1;
    gps << ", '" << j0 << "' " << j0;
    for (int j = j0+1; j <= j0+J; ++j) {
        gps << ", '" << j << "' " << j;
    }
    gps << ")" << std::endl;
    gps << "plot " << j0-1.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();



    /*
    T maxCoeffSca = -1.0;
    T maxCoeffWav = -1.0;
    int j = (*coeff.begin()).second.j, k = (*coeff.begin()).second.k;
    int j0 = j;
    int J  = j;
    T a_sca = 5000.0, a_wav = 5000.0;
    T b_sca = -5000.0, b_wav = -5000.0;
    
    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        j = (*it).second.j; k = (*it).second.k;
        j0 = std::min(j0, j);
        J  = std::max(J, j);

        if ((*it).second.xtype == XBSpline) {
            maxCoeffSca = std::max(maxCoeffSca, fabs((*it).first));
            a_sca = std::min(a_sca, basis.mra.phi.support(j,k).l1);
            b_sca = std::max(b_sca, basis.mra.phi.support(j,k).l2);
        }
        else {
            maxCoeffWav = std::max(maxCoeffWav, fabs((*it).first));
            a_wav = std::min(a_wav, basis.psi.support(j,k).l1);
            b_wav = std::max(b_wav, basis.psi.support(j,k).l2);
        }
    }
    
    T maxCoeff = std::max(maxCoeffWav,maxCoeffSca);
    T l1_sca = 0; // basis.mra.phi.support(0,0).l1;
    T l2_sca = 1; // basis.mra.phi.support(0,0).l2;
    
    T l1_wav = 0; // basis.psi.support(0,0).l1;
    T l2_wav = 1; // basis.psi.support(0,0).l2;

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        T lineWidth = 0.1;
        T ctr, fromX, toX, fromY, toY;
        T color = 0.0;

        if ((*it).second.xtype==XBSpline) {
            //int k1 = ceil(pow2i<T>((*it).second.j)*a_sca - l1_sca), k2 = floor(pow2i<T>((*it).second.j)*b_sca - l2_sca);
            //int N = k2 - k1 + 1;
            //fromX = a_sca + ((*it).second.k-k1)*(b_sca-a_sca)/(T)N;
            //toX   = a_sca + ((*it).second.k-k1+1)*(b_sca-a_sca)/(T)N;
            
            int j = (*it).second.j;
            int k = (*it).second.k;

            fromX = basis.mra.phi.support(j,k).l1;
            toX   = basis.mra.phi.support(j,k).l2;


            fromY = (*it).second.j-1.5;
            toY   = (*it).second.j-0.5;
            color = fabs((*it).first) / maxCoeffSca;
        } else {
            long int k1 = ceil(pow2i<T>((*it).second.j)*a_wav - l1_wav), k2 = floor(pow2i<T>((*it).second.j)*b_wav - l2_wav);
            long int N = k2 - k1 + 1;
            fromX = a_wav + ((*it).second.k-k1)*(b_wav-a_wav)/(T)N;
            toX   = fromX + std::max((b_wav-a_wav)/N,0.1);  //was 0.05

            fromY = (*it).second.j-0.5;
            toY   = (*it).second.j+0.5;
            color = fabs((*it).first) / maxCoeffWav;

        }

        gps << "set object rectangle from " << fromX << ", " << fromY
            << " to " << toX << "," << toY << " fc rgb";
        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";
        gps << " linewidth " << 0.1 << " fillstyle solid" << std::endl;

    }

    gps << "set xrange[0:1]" << std::endl;
    gps << "set yrange[" << j0-1.5 << ":" << J+0.5 << "]" << std::endl;
    //gps << "set xtics 1" << endl;
    gps << "set ytics ('" << j0 << "' " << j0-1;
    gps << ", '" << j0 << "' " << j0;
    for (int j = j0+1; j <= J; ++j) {
        gps << ", '" << j << "' " << j;
    }
    gps << ")" << std::endl;
    gps << "plot " << j0-1.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();
    */
}



template <typename T>
void
plotCoeff(const Coefficients<AbsoluteValue,T,Index1D > &coeff, const Basis<T,Primal,R,CDF> &basis, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator const_it;
    if (coeff.size() == 0) {
        return;
    }
    std::cout << "plotCoeff was called!" << std::endl;

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps enh color; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    T maxCoeffSca = -1.0;
    T maxCoeffWav = -1.0;
    int j = (*coeff.begin()).second.j, k = (*coeff.begin()).second.k;
    int j0 = j;
    int J  = j;
    T a_sca = 5000.0, a_wav = 5000.0;
    T b_sca = -5000.0, b_wav = -5000.0;
    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        j = (*it).second.j; k = (*it).second.k;
        j0 = std::min(j0, j);
        J  = std::max(J, j);
        if ((*it).second.xtype == XBSpline) {
            maxCoeffSca = std::max(maxCoeffSca, fabs((*it).first));
            a_sca = std::min(a_sca, basis.mra.phi.support(j,k).l1);
            b_sca = std::max(b_sca, basis.mra.phi.support(j,k).l2);
        }
        else {
            maxCoeffWav = std::max(maxCoeffWav, fabs((*it).first));
            a_wav = std::min(a_wav, basis.psi.support(j,k).l1);
            b_wav = std::max(b_wav, basis.psi.support(j,k).l2);
        }
    }
    T maxCoeff = std::max(maxCoeffWav,maxCoeffSca);
    T l1_sca = basis.mra.phi.support(0,0).l1, l2_sca = basis.mra.phi.support(0,0).l2;
    T l1_wav = basis.psi.support(0,0).l1,     l2_wav = basis.psi.support(0,0).l2;

    for (const_it it = coeff.begin(); it != coeff.end(); ++it) {
        T lineWidth = 0.1;
        T ctr, fromX, toX, fromY, toY;
        T color = 0.0;

        if ((*it).second.xtype==XBSpline) {
            int k1 = ceil(pow2i<T>((*it).second.j)*a_sca - l1_sca), k2 = floor(pow2i<T>((*it).second.j)*b_sca - l2_sca);
            int N = k2 - k1 + 1;
            fromX = a_sca + ((*it).second.k-k1)*(b_sca-a_sca)/(T)N;
            toX   = a_sca + ((*it).second.k-k1+1)*(b_sca-a_sca)/(T)N;

            fromY = (*it).second.j-1.5;
            toY   = (*it).second.j-0.5;
            color = fabs((*it).first) / maxCoeffSca;
        }

        else {
            long int k1 = ceil(pow2i<T>((*it).second.j)*a_wav - l1_wav), k2 = floor(pow2i<T>((*it).second.j)*b_wav - l2_wav);
            long int N = k2 - k1 + 1;
            fromX = a_wav + ((*it).second.k-k1)*(b_wav-a_wav)/(T)N;
            toX   = fromX + std::max((b_wav-a_wav)/N,0.1);  //was 0.05

            fromY = (*it).second.j-0.5;
            toY   = (*it).second.j+0.5;
            color = fabs((*it).first) / maxCoeffWav;

        }

        gps << "set object rectangle from " << fromX << ", " << fromY
            << " to " << toX << "," << toY << " fc rgb";
        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";
        gps << " linewidth " << 0.1 << " fillstyle solid" << std::endl;

    }

    gps << "set xrange["<< std::min(a_sca, a_wav) <<":"<< std::max(b_sca, b_wav) <<"]" << std::endl;
    gps << "set yrange[" << j0-1.5 << ":" << J+0.5 << "]" << std::endl;
    //gps << "set xtics 1" << endl;
    gps << "set ytics ('" << j0 << "' " << j0-1;
    gps << ", '" << j0 << "' " << j0;
    for (int j = j0+1; j <= J; ++j) {
        gps << ", '" << j << "' " << j;
    }
    gps << ")" << std::endl;
    gps << "plot " << j0-1.5 << " with lines linecolor rgb 'black' notitle" << std::endl;
    gps << "reset; set terminal pop" << std::endl;
    gps.close();
}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x, const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_coeff_abs_it;

    std::stringstream gpsFilename;
    gpsFilename << filename << ".gps";
    std::ofstream gps(gpsFilename.str().c_str());

    T h = 5e-2;

    gps << "reset" << std::endl;
    gps << "set terminal postscript eps color enh; set output '" << filename << ".eps'" << std::endl;
    gps << "set palette color; set colorbox vertical" << std::endl;

    const_coeff_abs_it first_element = coeff.begin();
    T max_value = fabs( (*first_element).first );
    T min_x=10000., max_x=-10000., min_y=10000., max_y=-10000.;
    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, k1=(*it).second.index1.k, j2=(*it).second.index2.j, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;
        
        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        min_x = std::min(min_x,x); max_x = std::max(max_x,x);
        min_y = std::min(min_y,y); max_y = std::max(max_y,y);
    }

    T ratio = (max_y-min_y)/(max_x-min_x);

    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, k1=(*it).second.index1.k, j2=(*it).second.index2.j, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;
     
        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        min_x = std::min(min_x,x); max_x = std::max(max_x,x);
        min_y = std::min(min_y,y); max_y = std::max(max_y,y);
        T color = fabs((*it).first)/max_value;

        if (ratio<1) {
            gps << "set object rectangle from " << x-h << ", " << y-h/ratio << " to " << x+h << "," << y+h/ratio << "fc rgb ";
        }
        else {
            gps << "set object rectangle from " << x-h/ratio << ", " << y-h << " to " << x+h/ratio << "," << y+h << "fc rgb ";
        }

        if (color > 0.5) gps << " 'black' ";
        else if ((0.5 >= color) && (color > 0.25)) gps << " 'purple' ";
        else if ((0.25 >= color) && (color > 0.125)) gps << " 'magenta' ";
        else if ((0.125 >= color) && (color > 0.0625)) gps << " 'red' ";
        else if ((0.0625 >= color) && (color > 0.03125)) gps << " 'orangered' ";
        else if ((0.03125 >= color) && (color > 0.015625)) gps << " 'orange' ";
        else if ((0.015625 >= color) && (color > 0.0078125)) gps << " 'yellow' ";
        else gps << " 'grey' ";

        gps << " fs solid 1.0" << std::endl;
    }
    gps << "set xlabel 'x'" << std::endl;
    gps << "set ylabel 'y'" << std::endl;
    gps << "set xrange[" << min_x-h << ":" << max_x+h << "]" << std::endl;
    gps << "set yrange[" << min_y-h << ":" << max_y+h << "]" << std::endl;
    gps << "plot " << std::min(min_x-h,min_y-h) << " w l lc rgb 'black' notitle " << std::endl;
    gps.close();
}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<AbsoluteValue,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_coeff_abs_it;


    std::stringstream dataFilename;
    dataFilename << filename << ".dat";
    std::ofstream data(dataFilename.str().c_str());

    for (const_coeff_abs_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).second.index1.j, k1=(*it).second.index1.k, j2=(*it).second.index2.j, k2=(*it).second.index2.k;
        XType type1=(*it).second.index1.xtype, type2=(*it).second.index2.xtype;

        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        data << x << " " << y << std::endl;
    }
    data.close();

}

template <typename T, typename Index, typename Basis_x, typename Basis_y>
void
plotScatterCoeff2D(const Coefficients<Lexicographical,T,Index> &coeff, const Basis_x &basis_x,
                   const Basis_y &basis_y, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;


    std::stringstream dataFilename;
    dataFilename << filename << ".dat";
    std::ofstream data(dataFilename.str().c_str());

    for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
        int j1=(*it).first.index1.j, k1=(*it).first.index1.k, j2=(*it).first.index2.j, k2=(*it).first.index2.k;
        XType type1=(*it).first.index1.xtype, type2=(*it).first.index2.xtype;

        //center of the support
        double x = 0.5*(basis_x.generator(type1).support(j1,k1).l2 + basis_x.generator(type1).support(j1,k1).l1);
        double y = 0.5*(basis_y.generator(type2).support(j2,k2).l2 + basis_y.generator(type2).support(j2,k2).l1);

        data << x << " " << y << std::endl;
    }
    data.close();

}


}  // namespace lawa

