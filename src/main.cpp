
#include <iostream>
#include <gtest/gtest.h>
#include "hw1.h"

int main(int argc, char **argv)
{
    if (true) // make false to run unit-tests
    {
        Matrix m{{6,1,1},{4,-2,5},{2,8,7}};
        // for(int i=1;i<=3;i++){
        //     for(int j=1;j<=3;j++){
        //         m[i-1][j-1]=i*j;
        //     }
        // }
        // Matrix c=algebra::minorr(m,1,1);
        // for(int i=0;i<2;i++){
        //     for(int j=0;j<2;j++){
        //         std::cout<<c[i][j]<<" ";
        //     }
        //     std::cout<<std::endl;
        // }
        std::cout<<algebra::determinant(m);
    }
    else
    {
        ::testing::InitGoogleTest(&argc, argv);
        std::cout << "RUNNING TESTS ..." << std::endl;
        int ret{RUN_ALL_TESTS()};
        if (!ret)
            std::cout << "<<<SUCCESS>>>" << std::endl;
        else
            std::cout << "FAILED" << std::endl;
    }
    return 0;
}