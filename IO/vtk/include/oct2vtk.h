//
// Created by milinda on 5/30/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief writes a given octree and variable array to vtk file.
 * The binary output is written in little endian.
*/
//

#ifndef SFCSORTBENCH_OCT2VTK_H
#define SFCSORTBENCH_OCT2VTK_H

#define FNAME_LENGTH 256
#define VTK_HEXAHEDRON 12

#define DENDRO_NODE_COORD_TYPE "UInt32"
#define DENDRO_NODE_ID_TYPE "UInt64"
#define DENDRO_NODE_VAR_INT "UInt32"
#define DENDRO_NODE_VAR_FLOAT "Float32"
#define DENDRO_NODE_VAR_DOUBLE "Float64"

#define DENDRO_FORMAT_ASCII "ascii"
#define DENDRO_FORMAT_BINARY "binary"




#include "TreeNode.h"
#include "mesh.h"
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "cencode.h"
#include "zlib.h"
#include <iostream>
#include <string>

namespace io
{
    namespace vtk
    {

        static FILE *fp = NULL;

        static void write_vtu_header(void)
        {
            fprintf(fp, "# vtk DataFile Version 2.0\n");
            fprintf(fp, "Dendro-5.0\n");

#ifdef DENDRO_VTU_ASCII
            fprintf(fp, "ASCII\n");
#else
            fprintf(fp, "BINARY\n");
#endif

        }


        static int vtk_write_compressed (FILE * vtkfile, char *numeric_data,size_t byte_length)
        {
            int                 retval, fseek1, fseek2;
            size_t              iz;
            size_t              blocksize, lastsize;
            size_t              theblock, numregularblocks, numfullblocks;
            size_t              header_entries, header_size;
            size_t              code_length, base_length;
            long                header_pos, final_pos;
            char               *comp_data, *base_data;
            uint32_t           *compression_header;
            uLongf              comp_length;
            base64_encodestate  encode_state;

            /* compute block sizes */
            blocksize = (size_t) (1 << 15);       /* 32768 */
            lastsize = byte_length % blocksize;
            numregularblocks = byte_length / blocksize;
            numfullblocks = numregularblocks + (lastsize > 0 ? 1 : 0);
            header_entries = 3 + numfullblocks;
            header_size = header_entries * sizeof (uint32_t);

            /* allocate compression and base64 arrays */
            code_length = 2 * std::max (blocksize, header_size) + 4 + 1;
            comp_data =(char*) malloc (code_length*sizeof(char));
            base_data =(char*) malloc(code_length*sizeof(char));

            /* figure out the size of the header and write a dummy */
            compression_header = (uint32_t*)malloc (header_entries*sizeof(uint32_t));
            compression_header[0] = (uint32_t) numfullblocks;
            compression_header[1] = (uint32_t) blocksize;
            compression_header[2] = (uint32_t)
                    (lastsize > 0 || byte_length == 0 ? lastsize : blocksize);
            for (iz = 3; iz < header_entries; ++iz) {
                compression_header[iz] = 0;
            }
            base64_init_encodestate (&encode_state);
            base_length = base64_encode_block ((char *) compression_header,
                                               header_size, base_data, &encode_state);
            base_length +=
                    base64_encode_blockend (base_data + base_length, &encode_state);
            assert (base_length < code_length);
            base_data[base_length] = '\0';
            header_pos = ftell (vtkfile);
            (void) fwrite (base_data, 1, base_length, vtkfile);

            /* write the regular data blocks */
            base64_init_encodestate (&encode_state);
            for (theblock = 0; theblock < numregularblocks; ++theblock) {
                comp_length = code_length;
                retval = compress2 ((Bytef *) comp_data, &comp_length,
                                    (const Bytef *) (numeric_data + theblock * blocksize),
                                    (uLong) blocksize, Z_BEST_COMPRESSION);
                assert (retval == Z_OK);
                compression_header[3 + theblock] = comp_length;
                base_length = base64_encode_block (comp_data, comp_length,
                                                   base_data, &encode_state);
                assert (base_length < code_length);
                base_data[base_length] = '\0';
                (void) fwrite (base_data, 1, base_length, vtkfile);
            }

            /* write odd-sized last block if necessary */
            if (lastsize > 0) {
                comp_length = code_length;
                retval = compress2 ((Bytef *) comp_data, &comp_length,
                                    (const Bytef *) (numeric_data + theblock * blocksize),
                                    (uLong) lastsize, Z_BEST_COMPRESSION);
                assert (retval == Z_OK);
                compression_header[3 + theblock] = comp_length;
                base_length = base64_encode_block (comp_data, comp_length,
                                                   base_data, &encode_state);
                assert (base_length < code_length);
                base_data[base_length] = '\0';
                (void) fwrite (base_data, 1, base_length, vtkfile);
            }

            /* write base64 end block */
            base_length = base64_encode_blockend (base_data, &encode_state);
            assert (base_length < code_length);
            base_data[base_length] = '\0';
            (void) fwrite (base_data, 1, base_length, vtkfile);

            /* seek back, write header block, seek forward */
            final_pos = ftell (vtkfile);
            base64_init_encodestate (&encode_state);
            base_length = base64_encode_block ((char *) compression_header,
                                               header_size, base_data, &encode_state);
            base_length +=
                    base64_encode_blockend (base_data + base_length, &encode_state);
            assert (base_length < code_length);
            base_data[base_length] = '\0';
            fseek1 = fseek (vtkfile, header_pos, SEEK_SET);
            (void) fwrite (base_data, 1, base_length, vtkfile);
            fseek2 = fseek (vtkfile, final_pos, SEEK_SET);

            /* clean up and return */
            free (compression_header);
            free (comp_data);
            free (base_data);
            if (fseek1 != 0 || fseek2 != 0 || ferror (vtkfile)) {
                return -1;
            }
            return 0;
        }


        static int vtk_write_binary (FILE * vtkfile, char *numeric_data, size_t byte_length)
        {

#ifdef DENDRO_VTU_ZLIB
            return vtk_write_compressed(vtkfile,numeric_data,byte_length);
#else
            size_t              chunks, chunksize, remaining, writenow;
            size_t              code_length, base_length;
            uint32_t            int_header;
            char                *base_data;
            base64_encodestate  encode_state;

            /* VTK format used 32bit header info */
            assert (byte_length <= (size_t) UINT32_MAX);

            /* This value may be changed although this is not tested with VTK */
            chunksize = (size_t) 1 << 15; /* 32768 */
            int_header = (uint32_t) byte_length;

            /* Allocate sufficient memory for base64 encoder */
            code_length = 2 * std::max (chunksize, sizeof (int_header));
            code_length = std::max (code_length, (size_t)4) + 1;
            base_data = (char*)calloc(code_length,sizeof(char));// (code_length*sizeof(char));

            base64_init_encodestate (&encode_state);
            base_length =base64_encode_block ((char *) &int_header, sizeof (int_header), base_data,&encode_state);
            assert (base_length < code_length);
            base_data[base_length] = '\0';
            (void) fwrite (base_data, 1, base_length, vtkfile);

            chunks = 0;
            remaining = byte_length;
            while (remaining > 0) {
                writenow = std::min (remaining, chunksize);
                base_length = base64_encode_block (numeric_data + chunks * chunksize,writenow, base_data, &encode_state);
                assert (base_length < code_length);
                base_data[base_length] = '\0';
                (void) fwrite (base_data, 1, base_length, vtkfile);
                remaining -= writenow;
                ++chunks;
            }

            base_length = base64_encode_blockend (base_data, &encode_state);
            assert (base_length < code_length);
            base_data[base_length] = '\0';
            (void) fwrite (base_data, 1, base_length, vtkfile);

            free(base_data);
            if (ferror (vtkfile)) {
                return -1;
            }
            return 0;
#endif


        }


        static std::string getFileName(const std::string& s) {

            char sep = '/';

#ifdef _WIN32
            sep = '\\';
#endif

            size_t i = s.rfind(sep, s.length());
            if (i != std::string::npos) {
                return(s.substr(i+1, s.length() - i));
            }

            return(s);
        }



        /**
         *@breif Writes the given mesh to a binary vtk (legacy format) file.
         * @param [in] pMesh: input mesh
         * @param [in] fPrefix: vtk file prefix
         * @param [in] numVars: Number of vertex point variables.
         * @param [in] varNames: list of variable names
         * @param [in] vars: double ** pointer to the varaibles.
         * */
        void mesh2vtk(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,double * filedData,unsigned int numPointdata, const char **pointDataNames, double **pointData);

        /**
        *@breif Writes the given mesh to a binary vtu (in xml format) file.
        * @param [in] pMesh: input mesh
        * @param [in] fPrefix: vtk file prefix
        * @param [in] numVars: Number of vertex point variables.
        * @param [in] varNames: list of variable names
        * @param [in] vars: double ** pointer to the varaibles.
        * */
        void mesh2vtu(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData);


        /**
        *@breif Writes the given octree to a binary vtu (in xml format) file.
        * @param [in] pNodes: input nodes
        * @param [in] fPrefix: vtk file prefix
        * @param [in] comm : mpi communicator.
        * */
        void oct2vtu(const ot::TreeNode *pNodes,const unsigned int nSize,const char* fPrefix, MPI_Comm comm);


        /**
        *@breif Writes the given mesh to a binary vtu (in xml format) file (Note that this will write only the coarse octree. ).
        * @param [in] pMesh: input mesh
        * @param [in] fPrefix: vtk file prefix
        * @param [in] numVars: Number of vertex point variables.
        * @param [in] varNames: list of variable names
        * @param [in] vars: double ** pointer to the varaibles.
        * */
        void mesh2vtuCoarse(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData);

        /**
        *@breif Writes the given mesh to a binary vtu (in xml format) file (Note that this will write only the coarse octree. ).
        * @param [in] pMesh: input mesh
        * @param [in] fPrefix: vtk file prefix
        * @param [in] numVars: Number of vertex point variables.
        * @param [in] varNames: list of variable names
        * @param [in] vars: double ** pointer to the varaibles.
        * */
        void mesh2vtuFine(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData);
        
        
        /**
         * @brief computes the wavelet coefficients for specified varIDs and write them and write them along with the mesh. 
         * @param [in] mesh: input mesh
         * @param [in] zippedVec zip version of the variable array. 
         * @param [in] unzippedVec: unzip version of the variable array. 
         * @param [in] varIds: variable IDs to compute the wavelets on
         * @param [in] numVars: number of variables
         * @param [in] fileName: output file name
         *
         */
        
        void waveletsToVTU(ot::Mesh* mesh,const double ** zippedVec, const double **unzippedVec,const unsigned int * varIds,const unsigned int numVars,const char* fileName);
        


    }
}



#endif //SFCSORTBENCH_OCT2VTK_H
