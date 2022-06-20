function [out] = grad_x(in,dim,dx)
    out=zeros(dim);
    idx=1.0/dx;
    idx_by_2=0.5*idx;
    idx_by_12=idx/12.0;
    
    for k=1:(dim(3))
        for j=1:(dim(2))
            
            for i=3: (dim(1)-2)
                out(i,j,k)= (in(i-2,j,k) -8.0*in(i-1,j,k)+8.0*in(i+1,j,k)-in(i+2,j,k))*idx_by_12;
            end
            
            pp=k*dim(2)*dim(1)+j*dim(1)+1;
            out(1,j,k) = (-3.0*in(1,j,k)+ 4.0*in(i+1,j,k)-in(i+2,j,k))*idx_by_2;
            
            pp=k*dim(2)*dim(1)+j*dim(1)+2;
            out(2,j,k) = (-in(1,j,k)+in(3,j,k))*idx_by_2;
            
            
            pp=k*dim(2)*dim(1)+j*dim(1)+dim(1)-2;
            out(dim(1)-2,j,k) = (-in(dim(1)-3,j,k) +in(dim(1)-1,j,k))*idx_by_2;
            
            pp=k*dim(2)*dim(1)+j*dim(1)+dim(1)-1;
            out(dim(1)-1,j,k) = (in(dim(1)-2,j,k) -4.0*in(dim(1)-1,j,k) + 3.0*in(dim(1)-1,j,k))*idx_by_2;
                        
        end
    end
    
    



end