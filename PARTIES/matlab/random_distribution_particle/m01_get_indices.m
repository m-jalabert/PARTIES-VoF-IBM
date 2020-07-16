%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% get indices %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i,iii,j,jjj,k,kkk] = get_indices(PX,iii,PY,jjj,PZ,kkk)

    i = 0;
    j = 0;
    k = 0;    
    if iii == 0
      i   = i-1;
      iii = PX;
    end  
    if iii == PX+1
      i   = i+1;            
      iii = 1;
    end   
    if jjj == 0
      j   = j-1;            
      jjj = PY;
    end   
    if jjj == PY+1
      j   = j+1;            
      jjj = 1;
    end     
    if kkk == 0
      k   = k-1;           
      kkk = PZ;
    end    
    if kkk == PZ+1
      k   = k+1;          
      kkk = 1;
    end 
end