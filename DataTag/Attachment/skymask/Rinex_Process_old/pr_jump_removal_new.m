function pr1_R = pr_jump_removal_new(pr1_R, threshold)


pr1_R0 = nan(size(pr1_R,1),size(pr1_R,2));
for idt = 1:size(pr1_R,2)
    for prn = 1:size(pr1_R,1)
        if pr1_R(prn,idt)~= 0 
            pr1_R0(prn,idt) = pr1_R(prn,idt);
        end
    end
end

pr_ed = nan(size(pr1_R,1),size(pr1_R,2));
for idt = 2:size(pr1_R,2)
    for prn = 1:size(pr1_R,1)
        if pr1_R(prn,idt)~=0 %&& pr1_R(prn,idt-1)~=0
            pr_ed(prn,idt) = pr1_R(prn,idt) - pr1_R(prn,idt-1);
        end
    end
end  

% flag_pr = nan(size(pr1_R,1),size(pr1_R,2));
% for prn = 1:size(pr1_R,1)
%     if pr1_R(prn,1)~=0
%         temp_flag = 1;
%     else
%         temp_flag = 0;
%     end
%     flag_pr(prn,1) = temp_flag;
%     for idt = 2:size(pr1_R,2)
%         if abs(pr_ed(prn,idt)-median(pr_ed(prn,:),'omitnan'))>10000000
%             temp_flag = 1;
%         elseif abs(pr_ed(prn,idt)-median(pr_ed(prn,:),'omitnan'))>threshold && idt > 3 && pr_ed(prn,idt-2)==1
%             temp_flag = 2;
%         elseif isnan(pr_ed(prn,idt))
%             temp_flag = 0;
%         end
%         flag_pr(prn,idt) = temp_flag;
%     end
% end
% 
% for prn = 101:size(pr1_R,1)
%     for idt = 1:size(pr1_R,2)
%         if flag_pr(prn,idt)==1
%             pr1_R(prn,idt) = 0;
%         end
%     end
% end


flag_pr = nan(size(pr1_R,1),size(pr1_R,2));
for prn = 1:size(pr1_R,1)
    if pr1_R(prn,1)~=0
        temp_flag = 1;
    else
        temp_flag = 0;
    end
    flag_pr(prn,1) = temp_flag;
    for idt = 2:size(pr1_R,2)
        if isnan(pr_ed(prn,idt)) || pr1_R(prn,idt)>40000000
            temp_flag = 0;
        elseif abs(pr_ed(prn,idt)-median(pr_ed(prn,:),'omitnan'))>10000000
            temp_flag = 1;
        elseif abs(pr_ed(prn,idt)-median(pr_ed(prn,:),'omitnan'))>threshold && idt > 3 && sum(flag_pr(prn,idt-2:idt-1))==2
            temp_flag = 2;
        end
        flag_pr(prn,idt) = temp_flag;
    end
end

for prn = 101:size(pr1_R,1)
    for idt = 1:size(pr1_R,2)
        if flag_pr(prn,idt)==1 || flag_pr(prn,idt)==0
            pr1_R(prn,idt) = 0;
        end
    end
end