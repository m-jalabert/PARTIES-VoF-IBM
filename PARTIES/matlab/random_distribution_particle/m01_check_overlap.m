%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% check overlap %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fail,p] = check_overlap(xp,yp,zp,rp,x,y,z,r,p)
    
    % particle - particle
    fail = 0;
    for q = 1: p;
       if p~=q
          dis = sqrt((xp-x(q))^2+(yp-y(q))^2+(zp-z(q))^2);           
%               fprintf('xp : %g yp : %g zp : %g %d  \n' , xp,   yp,   zp,   dis) 
%               fprintf('xp : %g yp : %g zp : %g %d\n\n' , x(q), y(q), z(q), fail)            
           if (dis < rp + r(q));
              fail = 1;   
           end
       end      
    end
    
    if fail == 1
      p = p - 1;
    end
    
  

    

end