function [orbit,reach_idx,delete,reach] = passplan(rvr,nR,nt)
if rvr == 5
    addpath('Kushiyara_passplan')
    orbit_161=table2array(imp_orbitfile('passplan_t0161.txt'));
    orbit_230=table2array(imp_orbitfile('passplan_t0230.txt'));
    orbit_467=table2array(imp_orbitfile('passplan_t0467.txt'));
    orbit_536=table2array(imp_orbitfile('passplan_t0536.txt'));
    
    % 536:1
    % 467:1,2,3,4,5
    % 230:1,2,3,4,5
    % 161:4,5
    % reach 1: orbit 230,467,536
    % reach 2: orbit 230,467
    % reach 3: orbit 230,467
    % reach 4: orbit 161,230,467
    % reach 5: orbit 161,230,467
    
    orbit{1} = ceil(orbit_161);
    orbit{2} = ceil(orbit_230);
    orbit{3} = ceil(orbit_467);
    orbit{4} = ceil(orbit_536);
    reach_idx{1} = [2,3,4];
    reach_idx{2} = [2,3];
    reach_idx{3} = [2,3];
    reach_idx{4} = [1,2,3];
    reach_idx{5} = [1,2,3];
    
    for i = 1:nR
        tmp = [];
        for j = 1:length(reach_idx{i})
            tmp = [tmp;orbit{reach_idx{i}(j)}];
        end
        reach{i} = tmp;
        delete{i} = setdiff(1:nt,reach{i});
    end
    
    
else if rvr ==1
        addpath('Ari_passplan')
        orbit_189=table2array(imp_orbitfile('passplan_t0189.txt'));
        orbit_258=table2array(imp_orbitfile('passplan_t0258.txt'));
        orbit_536=table2array(imp_orbitfile('passplan_t0536.txt'));
        
        % 189:1-10
        % 258:1-10
        % 536:9,10
        orbit{1} = ceil(orbit_189);
        orbit{2} = ceil(orbit_258);
        orbit{3} = ceil(orbit_536);
        orbit_idx{1}  = 1:10;
        orbit_idx{2}  = 1:10;
        orbit_idx{3}  = 1:10;
        for j = 1:nR
            tmp  =  [];
            for i = 1:length(orbit_idx)
                idx = (find(orbit_idx{i}==j));
                if ~isempty(idx)
                    tmp = [tmp,i];
                end
            end
            reach_idx{j}=tmp;
        end
        for i = 1:nR
            tmp = [];
            for j = 1:length(reach_idx{i})
                tmp = [tmp;orbit{reach_idx{i}(j)}];
            end
            reach{i} = tmp;
            delete{i} = setdiff(1:nt,reach{i});
        end
    else if rvr == 2
            addpath('Bra_passplan')
            orbit_189=table2array(imp_orbitfile('passplan_t0189.txt'));
            orbit_230=table2array(imp_orbitfile('passplan_t0230.txt'));
            orbit_258=table2array(imp_orbitfile('passplan_t0258.txt'));
            orbit_467=table2array(imp_orbitfile('passplan_t0467.txt'));
            orbit_495=table2array(imp_orbitfile('passplan_t0495.txt'));
            orbit_536=table2array(imp_orbitfile('passplan_t0536.txt'));
            
            
            orbit{1} = ceil(orbit_189);
            orbit{2} = ceil(orbit_230);
            orbit{3} = ceil(orbit_258);
            orbit{4} = ceil(orbit_467);
            orbit{5} = ceil(orbit_495);
            orbit{6} = ceil(orbit_536);
            orbit_idx{1}  = 1:5;
            orbit_idx{2}  = 3:5;
            orbit_idx{3}  = 1;
            orbit_idx{4}  = 4:5;
            orbit_idx{5}  = 1:2;
            orbit_idx{6}  = 1:4;
            for j = 1:nR
                tmp  =  [];
                for i = 1:length(orbit_idx)
                    idx = (find(orbit_idx{i}==j));
                    if ~isempty(idx)
                        tmp = [tmp,i];
                    end
                end
                reach_idx{j}=tmp;
            end
            for i = 1:nR
                tmp = [];
                for j = 1:length(reach_idx{i})
                    tmp = [tmp;orbit{reach_idx{i}(j)}];
                end
                reach{i} = tmp;
                delete{i} = setdiff(1:nt,reach{i});
            end
        else if rvr == 3
                addpath('Iowa_passplan')
                orbit_216=table2array(imp_orbitfile('passplan_t0216.txt'));
                orbit_259=table2array(imp_orbitfile('passplan_t0259.txt'));
                orbit_522=table2array(imp_orbitfile('passplan_t0522.txt'));
                
                
                orbit{1} = ceil(orbit_216);
                orbit{2} = ceil(orbit_259);
                orbit{3} = ceil(orbit_522);
                orbit_idx{1}  = 1:10;
                orbit_idx{2}  = 1:10;
                orbit_idx{3}  = 4:10;
                for j = 1:nR
                    tmp  =  [];
                    for i = 1:length(orbit_idx)
                        idx = (find(orbit_idx{i}==j));
                        if ~isempty(idx)
                            tmp = [tmp,i];
                        end
                    end
                    reach_idx{j}=tmp;
                end
                for i = 1:nR
                    tmp = [];
                    for j = 1:length(reach_idx{i})
                        tmp = [tmp;orbit{reach_idx{i}(j)}];
                    end
                    reach{i} = tmp;
                    delete{i} = setdiff(1:nt,reach{i});
                end
            else if rvr ==4
                    addpath('Jamu_passplan')
                    orbit_189=table2array(imp_orbitfile('passplan_t0189.txt'));
                    orbit_258=table2array(imp_orbitfile('passplan_t0258.txt'));
                    orbit_495=table2array(imp_orbitfile('passplan_t0495.txt'));
                    orbit_536=table2array(imp_orbitfile('passplan_t0536.txt'));
                    
                    
                    orbit{1} = ceil(orbit_189);
                    orbit{2} = ceil(orbit_258);
                    orbit{3} = ceil(orbit_495);
                    orbit{4} = ceil(orbit_536);
                    orbit_idx{1}  = 1:4;
                    orbit_idx{2}  = 1:4;
                    orbit_idx{3}  = 1:4;
                    orbit_idx{4}  = 3:4;
                    for j = 1:nR
                        tmp  =  [];
                        for i = 1:length(orbit_idx)
                            idx = (find(orbit_idx{i}==j));
                            if ~isempty(idx)
                                tmp = [tmp,i];
                            end
                        end
                        reach_idx{j}=tmp;
                    end
                    for i = 1:nR
                        tmp = [];
                        for j = 1:length(reach_idx{i})
                            tmp = [tmp;orbit{reach_idx{i}(j)}];
                        end
                        reach{i} = tmp;
                        delete{i} = setdiff(1:nt,reach{i});
                    end
                else if rvr ==6
                        addpath('MissDown_passplan')
                        orbit_009=table2array(imp_orbitfile('passplan_t0009.txt'));
                        orbit_272=table2array(imp_orbitfile('passplan_t0272.txt'));
                        idx1 = find(diff(orbit_009)<0)+1;
                        idx2 = find(diff(orbit_272)<0)+1;
                        orbit_009(idx1:end) = orbit_009(idx1:end)+365;
                        orbit_272(idx2:end) = orbit_272(idx2:end)+365;
                        
                        orbit{1} = ceil(orbit_009);
                        orbit{2} = ceil(orbit_272);
                        orbit_idx{1}  = 1:4;
                        orbit_idx{2}  = 1:5;
                        for j = 1:nR
                            tmp  =  [];
                            for i = 1:length(orbit_idx)
                                idx = (find(orbit_idx{i}==j));
                                if ~isempty(idx)
                                    tmp = [tmp,i];
                                end
                            end
                            reach_idx{j}=tmp;
                        end
                        for i = 1:nR
                            tmp = [];
                            for j = 1:length(reach_idx{i})
                                tmp = [tmp;orbit{reach_idx{i}(j)}];
                            end
                            reach{i} = tmp;
                            delete{i} = setdiff(1:nt,reach{i});
                        end
                    else if rvr ==7
                            addpath('MissUp_passplan')
                            orbit_009=table2array(imp_orbitfile('passplan_t0009.txt'));
                            idx1 = find(diff(orbit_009)<0)+1;
                            orbit_009(idx1:end) = orbit_009(idx1:end)+365;
                            
                            orbit{1} = ceil(orbit_009);
                            orbit_idx{1}  = 1:4;
                            for j = 1:nR
                                tmp  =  [];
                                for i = 1:length(orbit_idx)
                                    idx = (find(orbit_idx{i}==j));
                                    if ~isempty(idx)
                                        tmp = [tmp,i];
                                    end
                                end
                                reach_idx{j}=tmp;
                            end
                            for i = 1:nR
                                tmp = [];
                                for j = 1:length(reach_idx{i})
                                    tmp = [tmp;orbit{reach_idx{i}(j)}];
                                end
                                reach{i} = tmp;
                                delete{i} = setdiff(1:nt,reach{i});
                            end
                        else if rvr ==8
                                addpath('Ohio1_passplan')
                                orbit_160=table2array(imp_orbitfile('passplan_t0160.txt'));
                                orbit_175=table2array(imp_orbitfile('passplan_t0175.txt'));
                                orbit_438=table2array(imp_orbitfile('passplan_t0438.txt'));
                                orbit_453=table2array(imp_orbitfile('passplan_t0453.txt'));
                                
                                orbit{1} = floor(orbit_160);
                                orbit{2} = floor(orbit_175);
                                orbit{3} = floor(orbit_438);
                                orbit{4} = floor(orbit_453);
                                orbit_idx{1}  = 1:8;
                                orbit_idx{2}  = 1:10;
                                orbit_idx{3}  = 1:10;
                                orbit_idx{4}  = 8:10;
                                for j = 1:nR
                                    tmp  =  [];
                                    for i = 1:length(orbit_idx)
                                        idx = (find(orbit_idx{i}==j));
                                        if ~isempty(idx)
                                            tmp = [tmp,i];
                                        end
                                    end
                                    reach_idx{j}=tmp;
                                end
                                for i = 1:nR
                                    tmp = [];
                                    for j = 1:length(reach_idx{i})
                                        tmp = [tmp;orbit{reach_idx{i}(j)}];
                                    end
                                    reach{i} = tmp;
                                    delete{i} = setdiff(1:nt,reach{i});
                                end
                            else if rvr ==9
                                    addpath('Ohio2_passplan')
                                    orbit_160=table2array(imp_orbitfile('passplan_t0160.txt'));
                                    orbit_175=table2array(imp_orbitfile('passplan_t0175.txt'));
                                    orbit_438=table2array(imp_orbitfile('passplan_t0438.txt'));
                                    orbit_481=table2array(imp_orbitfile('passplan_t0481.txt'));
                                    
                                    orbit{1} = floor(orbit_160);
                                    orbit{2} = floor(orbit_175);
                                    orbit{3} = floor(orbit_438);
                                    orbit{4} = floor(orbit_481);
                                    orbit_idx{1}  = [1,2,4,5,6,7,8];
                                    orbit_idx{2}  = 1:7;
                                    orbit_idx{3}  = 8;
                                    orbit_idx{4}  = 1:3;
                                    for j = 1:nR
                                        tmp  =  [];
                                        for i = 1:length(orbit_idx)
                                            idx = (find(orbit_idx{i}==j));
                                            if ~isempty(idx)
                                                tmp = [tmp,i];
                                            end
                                        end
                                        reach_idx{j}=tmp;
                                    end
                                    for i = 1:nR
                                        tmp = [];
                                        for j = 1:length(reach_idx{i})
                                            tmp = [tmp;orbit{reach_idx{i}(j)}];
                                        end
                                        reach{i} = tmp;
                                        delete{i} = setdiff(1:nt,reach{i});
                                    end
                                else if rvr == 10
                                        addpath('Ohio3_passplan')
                                        orbit_160=table2array(imp_orbitfile('passplan_t0160.txt'));
                                        orbit_175=table2array(imp_orbitfile('passplan_t0175.txt'));
                                        orbit_466=table2array(imp_orbitfile('passplan_t0466.txt'));
                                        orbit_481=table2array(imp_orbitfile('passplan_t0481.txt'));
                                        
                                        
                                        orbit{1} = floor(orbit_160);
                                        orbit{2} = floor(orbit_175);
                                        orbit{3} = floor(orbit_466);
                                        orbit{4} = floor(orbit_481);
                                        orbit_idx{1}  = 6:14;
                                        orbit_idx{2}  = 4:14;
                                        orbit_idx{3}  = 1:14;
                                        orbit_idx{4}  = 1:14;
                                        for j = 1:nR
                                            tmp  =  [];
                                            for i = 1:length(orbit_idx)
                                                idx = (find(orbit_idx{i}==j));
                                                if ~isempty(idx)
                                                    tmp = [tmp,i];
                                                end
                                            end
                                            reach_idx{j}=tmp;
                                        end
                                        for i = 1:nR
                                            tmp = [];
                                            for j = 1:length(reach_idx{i})
                                                tmp = [tmp;orbit{reach_idx{i}(j)}];
                                            end
                                            reach{i} = tmp;
                                            delete{i} = setdiff(1:nt,reach{i});
                                        end
                                    else if rvr == 11
                                            addpath('Ohio4_passplan')
                                            orbit_466=table2array(imp_orbitfile('passplan_t0466.txt'));
                                            orbit_481=table2array(imp_orbitfile('passplan_t0481.txt'));
                                            
                                            
                                            orbit{1} = ceil(orbit_466);
                                            orbit{2} = ceil(orbit_481);
                                            orbit_idx{1}  = 1:5;
                                            orbit_idx{2}  = 1:4;
                                            for j = 1:nR
                                                tmp  =  [];
                                                for i = 1:length(orbit_idx)
                                                    idx = (find(orbit_idx{i}==j));
                                                    if ~isempty(idx)
                                                        tmp = [tmp,i];
                                                    end
                                                end
                                                reach_idx{j}=tmp;
                                            end
                                            for i = 1:nR
                                                tmp = [];
                                                for j = 1:length(reach_idx{i})
                                                    tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                end
                                                reach{i} = tmp;
                                                delete{i} = setdiff(1:nt,reach{i});
                                            end
                                        else if rvr == 12
                                                addpath('Ohio5_passplan')
                                                orbit_188=table2array(imp_orbitfile('passplan_t0188.txt'));
                                                orbit_203=table2array(imp_orbitfile('passplan_t0203.txt'));
                                                orbit_466=table2array(imp_orbitfile('passplan_t0466.txt'));
                                                orbit_481=table2array(imp_orbitfile('passplan_t0481.txt'));
                                                
                                                orbit{1} = ceil(orbit_188);
                                                orbit{2} = ceil(orbit_203);
                                                orbit{3} = ceil(orbit_466);
                                                orbit{4} = ceil(orbit_481);
                                                orbit_idx{1}  = 1:8;
                                                orbit_idx{2}  = 1:8;
                                                orbit_idx{3}  = 2:8;                                           orbit_idx{1}  = 1:5;
                                                orbit_idx{4}  = 1:8;
                                                for j = 1:nR
                                                    tmp  =  [];
                                                    for i = 1:length(orbit_idx)
                                                        idx = (find(orbit_idx{i}==j));
                                                        if ~isempty(idx)
                                                            tmp = [tmp,i];
                                                        end
                                                    end
                                                    reach_idx{j}=tmp;
                                                end
                                                for i = 1:nR
                                                    tmp = [];
                                                    for j = 1:length(reach_idx{i})
                                                        tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                    end
                                                    reach{i} = tmp;
                                                    delete{i} = setdiff(1:nt,reach{i});
                                                end
                                            else if rvr ==13
                                                    addpath('Ohio7_passplan')
                                                    orbit_188=table2array(imp_orbitfile('passplan_t0188.txt'));
                                                    orbit_203=table2array(imp_orbitfile('passplan_t0203.txt'));
                                                    orbit_494=table2array(imp_orbitfile('passplan_t0494.txt'));
                                                    orbit_509=table2array(imp_orbitfile('passplan_t0509.txt'));
                                                    
                                                    
                                                    orbit{1} = ceil(orbit_188);
                                                    orbit{2} = ceil(orbit_203);
                                                    orbit{3} = ceil(orbit_494);
                                                    orbit{4} = ceil(orbit_509);
                                                    orbit_idx{1}  = 7;
                                                    orbit_idx{2}  = 1:7;
                                                    orbit_idx{3}  = 1:7;
                                                    orbit_idx{4}  = 1:7;
                                                    for j = 1:nR
                                                        tmp  =  [];
                                                        for i = 1:length(orbit_idx)
                                                            idx = (find(orbit_idx{i}==j));
                                                            if ~isempty(idx)
                                                                tmp = [tmp,i];
                                                            end
                                                        end
                                                        reach_idx{j}=tmp;
                                                    end
                                                    for i = 1:nR
                                                        tmp = [];
                                                        for j = 1:length(reach_idx{i})
                                                            tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                        end
                                                        reach{i} = tmp;
                                                        delete{i} = setdiff(1:nt,reach{i});
                                                    end
                                                else if rvr ==14
                                                        addpath('Ohio8_passplan')
                                                        
                                                        orbit_203=table2array(imp_orbitfile('passplan_t0203.txt'));
                                                        orbit_216=table2array(imp_orbitfile('passplan_t0216.txt'));
                                                        orbit_494=table2array(imp_orbitfile('passplan_t0494.txt'));
                                                        orbit_509=table2array(imp_orbitfile('passplan_t0509.txt'));
                                                        
                                                        
                                                        orbit{1} = ceil(orbit_203);
                                                        orbit{2} = ceil(orbit_216);
                                                        orbit{3} = ceil(orbit_494);
                                                        orbit{4} = ceil(orbit_509);
                                                        orbit_idx{1}  = 5:8;
                                                        orbit_idx{2}  = 1;
                                                        orbit_idx{3}  = 1:8;
                                                        orbit_idx{4}  = 1:8;
                                                        for j = 1:nR
                                                            tmp  =  [];
                                                            for i = 1:length(orbit_idx)
                                                                idx = (find(orbit_idx{i}==j));
                                                                if ~isempty(idx)
                                                                    tmp = [tmp,i];
                                                                end
                                                            end
                                                            reach_idx{j}=tmp;
                                                        end
                                                        for i = 1:nR
                                                            tmp = [];
                                                            for j = 1:length(reach_idx{i})
                                                                tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                            end
                                                            reach{i} = tmp;
                                                            delete{i} = setdiff(1:nt,reach{i});
                                                        end
                                                    else if rvr ==15
                                                            addpath('Pada_passplan')
                                                            orbit_189=table2array(imp_orbitfile('passplan_t0189.txt'));
                                                            orbit_258=table2array(imp_orbitfile('passplan_t0258.txt'));
                                                            orbit_495=table2array(imp_orbitfile('passplan_t0495.txt'));
                                                            orbit_536=table2array(imp_orbitfile('passplan_t0536.txt'));
                                                            
                                                            
                                                            orbit{1} = ceil(orbit_189);
                                                            orbit{2} = ceil(orbit_258);
                                                            orbit{3} = ceil(orbit_495);
                                                            orbit{4} = ceil(orbit_536);
                                                            orbit_idx{1}  = 1:5;
                                                            orbit_idx{2}  = 1:5;
                                                            orbit_idx{3}  = 5;
                                                            orbit_idx{4}  = 2:5;
                                                            for j = 1:nR
                                                                tmp  =  [];
                                                                for i = 1:length(orbit_idx)
                                                                    idx = (find(orbit_idx{i}==j));
                                                                    if ~isempty(idx)
                                                                        tmp = [tmp,i];
                                                                    end
                                                                end
                                                                reach_idx{j}=tmp;
                                                            end
                                                            for i = 1:nR
                                                                tmp = [];
                                                                for j = 1:length(reach_idx{i})
                                                                    tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                                end
                                                                reach{i} = tmp;
                                                                delete{i} = setdiff(1:nt,reach{i});
                                                            end
                                                        else if rvr == 16
                                                                addpath('SeineDown_passplan')
                                                                orbit_113=table2array(imp_orbitfile('passplan_t0113.txt'));
                                                                orbit_292=table2array(imp_orbitfile('passplan_t0292.txt'));
                                                                orbit_391=table2array(imp_orbitfile('passplan_t0391.txt'));
                                                                
                                                                
                                                                orbit{1} = ceil(orbit_113);
                                                                orbit{2} = ceil(orbit_292);
                                                                orbit{3} = ceil(orbit_391);
                                                                orbit_idx{1}  = 1:6;
                                                                orbit_idx{2}  = [1,2,3,4,6];
                                                                orbit_idx{3}  = 1:6;
                                                                for j = 1:nR
                                                                    tmp  =  [];
                                                                    for i = 1:length(orbit_idx)
                                                                        idx = (find(orbit_idx{i}==j));
                                                                        if ~isempty(idx)
                                                                            tmp = [tmp,i];
                                                                        end
                                                                    end
                                                                    reach_idx{j}=tmp;
                                                                end
                                                                for i = 1:nR
                                                                    tmp = [];
                                                                    for j = 1:length(reach_idx{i})
                                                                        tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                                    end
                                                                    reach{i} = tmp;
                                                                    delete{i} = setdiff(1:nt,reach{i});
                                                                end
                                                                 else if rvr == 17
                                                                addpath('SeineUp_passplan')
                                                                orbit_085=table2array(imp_orbitfile('passplan_t0085.txt'));
                                                                orbit_292=table2array(imp_orbitfile('passplan_t0292.txt'));
                                                                orbit_391=table2array(imp_orbitfile('passplan_t0391.txt'));
                                                                orbit_570=table2array(imp_orbitfile('passplan_t0570.txt'));
                                                                
                                                                orbit{1} = ceil(orbit_085);
                                                                orbit{2} = ceil(orbit_292);
                                                                orbit{3} = ceil(orbit_391);
                                                                orbit{4} = ceil(orbit_570);
                                                                orbit_idx{1}  = 3:9;
                                                                orbit_idx{2}  = 1:9;
                                                                orbit_idx{3}  = 1:9;
                                                                orbit_idx{4}  = 3:9;
                                                                for j = 1:nR
                                                                    tmp  =  [];
                                                                    for i = 1:length(orbit_idx)
                                                                        idx = (find(orbit_idx{i}==j));
                                                                        if ~isempty(idx)
                                                                            tmp = [tmp,i];
                                                                        end
                                                                    end
                                                                    reach_idx{j}=tmp;
                                                                end
                                                                for i = 1:nR
                                                                    tmp = [];
                                                                    for j = 1:length(reach_idx{i})
                                                                        tmp = [tmp;orbit{reach_idx{i}(j)}];
                                                                    end
                                                                    reach{i} = tmp;
                                                                    delete{i} = setdiff(1:nt,reach{i});
                                                                end
                                                                     end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end