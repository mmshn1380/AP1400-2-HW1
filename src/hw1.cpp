#include "hw1.h"

namespace algebra{
    Matrix zeros(size_t n, size_t m){
        Matrix matrix(n , std::vector<double>(m , 0.0));
        return matrix;
    }
    Matrix ones(size_t n, size_t m){
        Matrix matrix(n , std::vector<double>(m , 1.0));
        return matrix;
    }
    Matrix random(size_t n, size_t m, double min, double max){
        Matrix matrix(n , std::vector<double>(m , 0.0));
        long int seed{std::chrono::system_clock::now().time_since_epoch().count()};
        std::minstd_rand0 generator (seed);
        for(size_t i{0} ; i<n ; i++)
            for(size_t j{0} ; j<m ; j++){
                matrix[i][j]=(static_cast<float>(generator())-static_cast<float>(generator.min()))/
                static_cast<float>(generator.max()-generator.min())*(max-min) + min;
            }
        return matrix;
    }
    Matrix multiply(const Matrix& matrix, double c){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx{matrix};
        for(size_t i{0} ; i<n ; i++)
            for(size_t j{0} ; j<m ; j++)
                matrixx[i][j]*=c;
        return matrixx;
    }
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
        size_t n1{matrix1.size()};
        size_t m1{matrix1[0].size()};
        size_t n2{matrix2.size()};
        size_t m2{matrix2[0].size()};
        if(m1!=n2){
            throw std::logic_error( "error in shape of matrix" );
        }
        Matrix mut(n1,std::vector<double>(m2,0.0));
        for(size_t i{0} ; i<n1 ; i++)
            for(size_t j{0} ; j<m2 ; j++){
                double sum{0};
                for(size_t k{0} ; j<m1 ; k++)
                    sum+=matrix1[i][k]+matrix2[k][j];
                mut[i][j]=sum;
            }
        
    }
    Matrix sum(const Matrix& matrix, double c){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx{matrix};
        for(size_t i{0} ; i<n ; i++)
            for(size_t j{0} ; j<m ; j++)
                matrixx[i][j]+=c;
        return matrixx;
    }
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2){
        size_t n1{matrix1.size()};
        size_t m1{matrix1[0].size()};
        size_t n2{matrix2.size()};
        size_t m2{matrix2[0].size()};
        if(n1!=n2 || m1!=m2){
            throw std::logic_error( "error in shape of matrix" );
        }
        Matrix matrix{matrix1};
        for(size_t i{0} ; i<n1 ; i++)
            for(size_t j{0} ; j<m1 ; j++)
                matrix[i][j]+=matrix2[i][j];
        return matrix;
        
    }
    Matrix transpose(const Matrix& matrix){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx(m , std::vector<double>(n , 0.0));
        for(size_t i{0} ; i<n ; i++)
            for(size_t j{0} ; j<m ; j++)
                matrixx[j][i]=matrix[i][j];
    }
    Matrix minorr(const Matrix& matrix, size_t n, size_t m){
        size_t n1{matrix.size()};
        size_t m1{matrix[0].size()};
        Matrix matrixx(n1-1 , std::vector<double>(m1-1 , 0.0));
        for(size_t i{0} ; i<n1 ; i++){
            for(size_t j{0} ; j<m1 ; j++){
                if(i!=n && j!=m){
                    matrixx[i<n?i:i-1][j<m?j:j-1]=matrix[i][j];
                }
            }
        }
        return matrixx;
    }
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){
        size_t n1{matrix1.size()};
        size_t m1{matrix1[0].size()};
        size_t n2{matrix2.size()};
        size_t m2{matrix2[0].size()};
        if(axis==0 && m1!=m2){
            throw std::logic_error( "error in shape of matrix" );
        }
        if(axis==1 && n1!=n2){
            throw std::logic_error( "error in shape of matrix" );
        }
        if(axis==0){
            Matrix matrixx(n1+n2 , std::vector<double>(m1 , 0.0));
            for(size_t i{0} ; i<n1 ; i++){
                for(size_t j{0} ; j<m1 ; j++){
                    matrixx[i][j]=matrix1[i][j];
                }
            }
            for(size_t i{0} ; i<n2 ; i++){
                for(size_t j{0} ; j<m2 ; j++){
                    matrixx[i+n1][j]=matrix2[i][j];
                }
            }
            return matrixx;
        }else{
            Matrix matrixx(n1 , std::vector<double>(m1+m2 , 0.0));
            for(size_t i{0} ; i<n1 ; i++){
                for(size_t j{0} ; j<m1 ; j++){
                    matrixx[i][j]=matrix1[i][j];
                }
            }
            for(size_t i{0} ; i<n2 ; i++){
                for(size_t j{0} ; j<m2 ; j++){
                    matrixx[i][j+m1]=matrix2[i][j];
                }
            }
            return matrixx;
        }
    }
    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx{matrix};
        for(size_t j{0} ; j<m ; j++){
            double save{matrixx[r1][j]};
            matrixx[r1][j]=matrixx[r2][j];
            matrixx[r2][j]=save;
        }
        return matrixx;
    }
    Matrix ero_multiply(const Matrix& matrix, size_t r, double c){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx{matrix};
        for(size_t j{0} ; j<m ; j++){
            matrixx[r][j]*=c;
        }
        return matrixx;
    }
    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        Matrix matrixx{matrix};
        for(size_t j{0} ; j<m ; j++){
            matrixx[r2][j]+=matrixx[r1][j]*c;
        }
        return matrixx;
    }
    Matrix inverse(const Matrix& matrix){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        if(n!=m){
            throw std::logic_error( "None square matrix does not have inverse" );
        }
        if(determinant(matrix)<0.00000001){
            throw std::logic_error( "Matrix with zero det have no inverse" );
        }
        Matrix adj{zeros(n,m)};
        for(size_t i{0} ; i<n ; i++){
            for(size_t j{0} ; j<m ; j++){
                adj[i][j]= ((i+j)%2==0 ? 1 : -1) * determinant(minorr(matrix,i,j));
            }
        }
        adj=transpose(adj);
        return multiply(adj,1/determinant(matrix));
    }
    Matrix upper_triangular(const Matrix& matrix){

    }
    double determinant(const Matrix& matrix){
        size_t n{matrix.size()};
        size_t m{matrix[0].size()};
        if(n!=m){
            throw std::logic_error( "None square matrix does not have determinal" );
        }
        if(n==1){
            return matrix[0][0];
        }
        double ans{};
        for(size_t j{};j<m;j++){
            ans+=(j%2==0 ? 1 : -1) * matrix[0][j] * determinant(minorr(matrix,0,j));
        }
        return ans;
    }
    void show(const Matrix& matrix){

    }
}