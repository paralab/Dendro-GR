//
// Created by milinda on 12/15/16.
//

#include "dendroIO.h"

namespace IO {
    int readPtsFromFile(char *filename, std::vector<double> &pts) {
        FILE *infile;
        size_t res;
        infile = fopen(filename, "rb");
        unsigned int temp;
        res = fread(&temp, sizeof(unsigned int), 1, infile);

        double *ptsTemp = NULL;

        if (temp) {
            ptsTemp = new double[3 * temp];
            assert(ptsTemp);
            res = fread(ptsTemp, sizeof(double), 3 * temp, infile);
        }

        fclose(infile);
        pts.resize(3 * temp);

        for (int i = 0; i < (3 * temp); i++) {
            pts[i] = ptsTemp[i];
        }//end for

        if (ptsTemp) {
            delete[] ptsTemp;
            ptsTemp = NULL;
        }

        // std::cout << __func__ << ": size " << temp << ", " << pts.size() << std::endl;
        return 1;
    }//end function

    int readDataPtsFromFile(char *filename, std::vector<double> &pts,
                            std::vector<double> &ptVals) {
        // file format:
        // 4 bytes (unsigned int?)  number of points N
        // 3*N*8 bytes coordinates of point (double);  X1 Y1 Z1 X2 Y2 Z2 ....
        // N*8 weights attached to points (double)
        //
        FILE *infile;
        size_t res;
        unsigned int temp;

        assert(sizeof(double) == 8);
        assert(sizeof(unsigned int) == 4);

        infile = fopen(filename, "rb");
        res = fread(&temp, sizeof(unsigned int), 1, infile);

        pts.resize(3 * temp);
        ptVals.resize(temp);

        res = fread(&(pts[0]), sizeof(double), 3 * temp, infile);
        res = fread(&(ptVals[0]), sizeof(double), temp, infile);

        fclose(infile);

        return 1;
    }//end function

    int writePtsToFile(char *filename, std::vector<double> &pts) {
        FILE *outfile = fopen(filename, "wb");
        unsigned int ptsLen = pts.size();
        double *ptsTemp = NULL;
        if (!pts.empty()) {
            ptsTemp = (&(*(pts.begin())));
        }

        if (ptsLen > 0) {
            unsigned int numPts = ptsLen / 3;
            fwrite(&numPts, sizeof(unsigned int), 1, outfile);
            fwrite(ptsTemp, sizeof(double), ptsLen, outfile);
        }
        fclose(outfile);
        return 1;
    }//end function


    int writeDataPtsToFile(char *filename, std::vector<double> &pts,
                           std::vector<double> &data) {
        FILE *outfile = fopen(filename, "wb");
        unsigned int ptsLen = pts.size();
        unsigned int numPts = ptsLen / 3;

        assert(sizeof(double) == 8);
        assert(sizeof(unsigned int) == 4);

        fwrite(&numPts, sizeof(unsigned int), 1, outfile);
        fwrite(&(pts[0]), sizeof(double), ptsLen, outfile);
        //fwrite(&(data[0]), sizeof(double),numPts,outfile);
        fclose(outfile);
        return 1;
    }//end function


    int readNodesFromFile(char *filename, std::vector<ot::TreeNode> &nodes) {
        int res;
        FILE *infile = fopen(filename, "r");
        unsigned int numNode;
        unsigned int dim, maxDepth;
        res = fscanf(infile, "%u", &dim);
        res = fscanf(infile, "%u", &maxDepth);
        res = fscanf(infile, "%u", &numNode);
        nodes.resize(numNode);

        for (unsigned int i = 0; i < nodes.size(); i++) {
            unsigned int x, y, z, d;
            res = fscanf(infile, "%u", &x);
            res = fscanf(infile, "%u", &y);
            res = fscanf(infile, "%u", &z);
            res = fscanf(infile, "%u", &d);
            nodes[i] = ot::TreeNode(x, y, z, d, dim, maxDepth);
        }
        fclose(infile);
        return 1;
    }//end function

    int writeNodesToFile(char *filename, const std::vector<ot::TreeNode> &nodes) {
        FILE *outfile = fopen(filename, "w");
        if (!nodes.empty()) {
            unsigned int dim = nodes[0].getDim();
            unsigned int maxDepth = nodes[0].getMaxDepth();
            fprintf(outfile, "%u %u\n", dim, maxDepth);
            fprintf(outfile, "%u\n", static_cast<unsigned int>(nodes.size()));
            for (unsigned int i = 0; i < nodes.size(); i++) {
                assert(nodes[i].getDim() == dim);
                assert(nodes[i].getMaxDepth() == maxDepth);
                fprintf(outfile, "%u %u %u %u\n", nodes[i].getX(), nodes[i].getY(), nodes[i].getZ(),
                        nodes[i].getLevel());
            }
        }
        fclose(outfile);
        return 1;
    }//end function
}