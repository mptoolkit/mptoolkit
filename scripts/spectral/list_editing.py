def get_next_smallest_w(finished_jobs_list,running_jobs_list,w0,configlist,k0):
    w_before='-1000'
    next_smaller=False
    for x in finished_jobs_list:
        x_split=x.split(',')
        if x_split[0]==k0:
            w_compare=x_split[1]
            if float(w_compare)<float(w0) and float(w_compare)>=float(w_before):
                next_smaller=x
                w_before=w_compare

    if next_smaller != False:
        for x in running_jobs_list:
            x_split=x.split(',')
            w_compare=x_split[1]
            w_next_smaller=(next_smaller.split(','))[1]
            if x_split[0]==k0:
                if float(w_compare)<float(w0) and float(w_compare)>float(w_next_smaller):
                    return False
        
    return next_smaller
        
    
def get_next_highest_w(finished_jobs_list,running_jobs_list,w0,configlist,k0):
    w_before='1000'
    next_higher=False
    for x in finished_jobs_list:
        x_split=x.split(',')
        if x_split[0]==k0:
            w_compare=x_split[1]
            if float(w_compare)>float(w0) and float(w_compare)<=float(w_before):
                next_higher=x
                w_before=w_compare
                
    if next_higher != False:
        for x in running_jobs_list:
            x_split=x.split(',')
            w_compare=x_split[1]
            w_next_higher=next_higher.split(',')[1]
            if x_split[0]==k0:
                if float(w_compare)>float(w0) and float(w_compare)<float(w_next_higher):
                    return False
        
    return next_higher
