import iovtk as iovtk
import filters as filters
import vtk as vtk
import argparse as argparse
from mpi4py import MPI

def main():
    comm = vtk.vtkMPIController()
    #c.SetGlobalController(None)
    rank = comm.GetLocalProcessId()
    npes = comm.GetNumberOfProcesses()
    #print ("rank %d of size %d" %(rank,npes))


    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--begin', help='begin step', required=True)
    parser.add_argument('-e','--end', help='end step', required=True)
    parser.add_argument('-f','--freq', help='io step frequency', required=True)
    parser.add_argument('-pvtu','--pvtu_prefix', help='pvtu prefix', required=True)
    #parser.add_argument('-img','--img_prefix', help='image prefix ', required=True)
    args = vars(parser.parse_args())
    args = parser.parse_args()



    step_begin=int(args.begin)
    step_end=int(args.end)
    step_freq=int(args.freq)
    file_prefix=args.pvtu_prefix
    #img_prefix=args.img_prefix

    '''step_begin=0
    step_end=20
    step_freq=10
    file_prefix='vtu/bssn_gr'
    img_prefix='img/img' 
    '''


    if(rank==0):
        print ("number of ranks: %d"%(npes))


    for step in range(step_begin,step_end,step_freq):
        fname = file_prefix + "_" + str(step) + ".pvtu"
        fslice = file_prefix + "slice_" + str(step) + ".pvtu"
        
        pvtuReader=iovtk.ReadPVTUFile(fname)
        gridSlice=filters.SliceFilter(pvtuReader)
        iovtk.WriteVTUFile(fslice,pvtuReader)

        print("slice written to %s " %fslice)
    
    return


if(__name__=="__main__"):
    main()

