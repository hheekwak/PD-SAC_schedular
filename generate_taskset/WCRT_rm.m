function [ R ] = WCRT_rm ( task )

    R = {}; flag = 1;
    for i = 1 : length(task)
        rq = 0;
        q = rq + 1;
        R(i).re(q) = task(i).C;
        while flag
            W_total = 0;
            for k = 1 : length(task)
                if k < i
                    W_total = W_total + ceil(R(i).re(q)/task(k).T)*task(k).C;
                end
            end
            
            if R(i).re(q) <= task(i).D && (task(i).C + W_total) > R(i).re(q)
                R(i).re(q) = task(i).C + W_total;
                rq = rq + 1;
                q = rq + 1;
                R(i).re(q) = R(i).re(q-1);
            elseif R(i).re(q) <= task(i).D && (task(i).C + W_total) <= R(i).re(q)
%                 R(i).re(q) = task(i).C + W_total;
                R(i).schedulable = 'schedulable';
                break;
            elseif R(i).re(q) > task(i).D
                R(i).schedulable = 'un-schedulable';
                break;
            end
        end
    end    

end

