do_umap <- function(V) {
    uwot::umap(
        X = V,
        n_threads = 6,
        n_neighbors = 30L,
        n_components = 2L,
        metric = 'cosine',
        n_epochs = NULL,
        learning_rate = 1.0,
        min_dist = 0.3,
        spread = 1.0,
        set_op_mix_ratio = 1.0,
        local_connectivity = 1L,
        repulsion_strength = 1,
        negative_sample_rate = 1,
        a = NULL,
        b = NULL,
        fast_sgd = FALSE,
        verbose = FALSE
    )    
    
}

do_scatter <- function (
    umap_use, meta_data, label_name, no_guides = FALSE, 
    do_labels = FALSE, nice_names, palette_use = tableau_color_pal()(10), 
    pt_size = 4, point_size = 0.5, pt_shape = '.',
    base_size = 12, do_points = TRUE, do_density = FALSE) {
    plt_df <- umap_use %>% data.frame() %>% cbind(meta_data) %>% 
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    if (!missing(nice_names)) {
        plt_df %<>% dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        plt_df[[label_name]] <- plt_df$nice_name
    }
    plt <- plt_df %>% ggplot(aes_string("X1", "X2", col = label_name, 
        fill = label_name)) + theme_tufte(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 16, size = 4)), alpha = FALSE) + 
        scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
        theme(plot.title = element_text(hjust = 0.5)) + labs(x = "UMAP 1", 
        y = "UMAP 2")
    if (do_points) 
        plt <- plt + geom_point(shape = pt_shape, size = point_size)
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    if (do_labels) {
        plt <- plt + geom_text_repel(data = data.table(plt_df)[, 
            .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
            label.size = NA, aes_string(label = label_name), 
            color = "black", size = pt_size, alpha = 1, segment.size = 0) + 
            guides(col = FALSE, fill = FALSE)
    }
    return(plt)
}

fig.size <- function (height, width) 
{
    options(repr.plot.height = height, repr.plot.width = width)
}

path_join <- function(...){
    c = paste(...,sep='/')
    return(c)      
}
makedirs <- function(path){
    dir.create(path,recursive = TRUE,mode = '0700')
}
removefile_or_path <- function(path_or_file){
    unlink(path_or_file,recursive = TRUE,force = TRUE)
}

mask_celltype <- function(meta_data, cell_type_use, mask_var_use="dataset", mask_var_values_use=c("jurkat", "t293")){
  cell_type = meta_data[[cell_type_use]]
  mask_var = meta_data[[mask_var_use]]
  mask_index = rep(FALSE,length(cell_type))
  logical_or = function(x) x[1]|x[2]
  for( mvu in mask_var_values_use){
    mask = mask_var == mvu
    mask_index = apply(cbind(mask_index,mask), 1,logical_or )  
  }
  fake_cell_type = as.character(cell_type)
  fake_cell_type[mask_index] = "unknown"
  if(is.factor(cell_type)){
    fake_cell_type = as.factor(fake_cell_type)
  }
  return(fake_cell_type)
}

get_embedding <- function(obj,pca_embedding){
    embedding <- t(as.matrix(obj$Z_corr))
    row.names(embedding ) <- row.names(pca_embedding)
    colnames(embedding ) <- colnames(pca_embedding )
    return(embedding)
}

lisi_mean_std = function(lisi_subtype_res){
    types = unique(lisi_subtype_res$type)
    keys = unique(lisi_subtype_res$key)
    #vals = matrix(0,nrow=length(types),ncol=length(keys))
    #counts = matrix(0,nrow=length(types),ncol=length(keys))
    means_stds =  matrix(0,nrow=length(types),ncol=length(keys)*2)
    rownames(means_stds ) = types
    means_stds_name = Reduce(c,lapply(keys, function(x) return(c(paste(x,"mean"), paste(x, 'std')))))
    colnames(means_stds ) = means_stds_name 
    i=1
    for(t in types){
        j=1
        for(k in keys){
            t_k_vector = lisi_subtype_res[["val"]][ lisi_subtype_res[['type']] ==t & lisi_subtype_res[['key']] == k]
            num_na = sum(is.na(t_k_vector))
            if(num_na>0){
                print(paste("here are ", num_na, "for ", t, "and",k))
            }
            means_stds[i,2*j-1] = mean(t_k_vector)
            means_stds[i,2*j] = sqrt(var(t_k_vector))
#             print(paste(t,k,vals[i,j]) )
            j=j+1
        }
        i=i+1
    }
    return( means_stds)
}

lisi_sum = function(lisi_subtype_res){
    types = unique(lisi_subtype_res$type)
    keys = unique(lisi_subtype_res$key)
    vals = matrix(0,nrow=length(types),ncol=length(keys))
    rownames(vals) = types
    colnames(vals) = keys
    i=1
    for(t in types){
        j=1
        for(k in keys){

            vals[i,j] = sum(lisi_subtype_res[["val"]][ lisi_subtype_res[['type']] ==t & lisi_subtype_res[['key']] == k])
#             print(paste(t,k,vals[i,j]) )
            j=j+1
        }
        i=i+1
    }
    return(vals)
}

onehot <- function(x) {
    data.frame(x) %>% 
        tibble::rowid_to_column("row_id") %>% 
        dplyr::mutate(dummy = 1) %>% 
        tidyr::spread(x, .data$dummy, fill = 0) %>% 
        dplyr::select(-.data$row_id) %>% as.matrix
}


# match score of R with the true Phi_C
get_fro_cluster_error_loss = function(obj,Phi_C)
{
    R = obj$R
    #Phi_C = harmonys2_obj$Phi_C
    M = R%*%t(Phi_C)
    K=dim(M)[1]
    C=dim(M)[2]

    Mp = diag(c(1./(M%*%matrix(1.,nrow=C,ncol=1)+1e-8))) %*% M 
    Psi_C = t(Mp) %*% R
    # to do normalize Psi_C
    diff = Psi_C - Phi_C
    error= norm(diff,'F')
    return(error)
}

summary_fro_cluster_error_loss <- function(objs,Phi_C){
    summary_table = data.frame(
        type=names(objs),
        fro_cluster_error_loss = Reduce(c, lapply(objs, function(obj) get_fro_cluster_error_loss(obj,Phi_C)                   
                                                )
                                       )
        )
    rownames( summary_table ) = NULL 
    return(summary_table)
}

# match score of R with the true Phi_C
get_entropy_error = function(obj,Phi_C)
{
    R = obj$R
    #Phi_C = harmonys2_obj$Phi_C
    M = R%*%t(Phi_C)
    K=dim(M)[1]
    C=dim(M)[2]

    Mp = diag(c(1./(M%*%matrix(1.,nrow=C,ncol=1)+1e-8))) %*% M 
    Psi_C = t(Mp) %*% R
    # to do normalize Psi_C
    kl = sum(Phi_C*log(Phi_C/(Psi_C+1e-8)),na.rm = T)
    return(kl)
}

summary_entropy_error <- function(objs,Phi_C){
    summary_table = data.frame(
        type=names(objs),
        fro_cluster_error_loss = Reduce(c, lapply(objs, function(obj) get_entropy_error(obj,Phi_C)                   
                                                )
                                       )
        )
    rownames( summary_table ) = NULL 
    return(summary_table)
}

                                      
                                
                                      
                                      
                                      

get_accuracy <- function(harmony_obj,Phi_C){
    R = harmony_obj$R
    #Phi_C = harmonys2_full_obj$Phi_C
    M = R%*%t(Phi_C)
    K=dim(M)[1]
    C=dim(M)[2]

    Mp = diag(c(1./( M%*%matrix(1.,nrow=C,ncol=1) + 1e-8 ) )) %*% M 
    Psi_C = t(Mp) %*% R
    batches_harmonys2 = apply(Psi_C,2, which.max)
    batches = apply(Phi_C,2, which.max)
    cave <- function(x, c1, c2) x[c1]== x[c2]
    X = cbind(batches,batches_harmonys2)
    colnames(X)<-c("x1","x2")
    res = apply(X, 1, cave,  c1 = "x1", c2 = c("x2"))
    accuracy = sum(res)*1./length(res)
    return(c(accuracy,sum(res),length(res)))
}

summary_accuracy <- function(objs,Phi_C){
    vals = Reduce(rbind, lapply(objs, function(obj) get_accuracy(obj,Phi_C)))
    colnames(vals) = c("accuracy","correct type", "total number of cells")
    summary_table = data.frame(
        type= names(objs)
        )
    summary_table = cbind(summary_table,vals)
    rownames( summary_table ) = NULL                           
    return(summary_table)
}

get_dincta_Phi_C_accuracy <- function(dincta_meta_data,cell_type="cell_type", 
                                      fake_cell_type = "fake_cell_type", cell_type_unknown = "unknown"){
    cell_type_real = as.character(dincta_meta_data[[cell_type]])
    cell_type_fake =   as.character(dincta_meta_data[[fake_cell_type]])
    unknow_index = cell_type_fake == cell_type_unknown
    cell_type_predict = dincta_meta_data[["cell_type_predict"]]
    cell_type_real_unknown = cell_type_real[unknow_index]
    cell_type_predict_unknown = cell_type_predict[unknow_index]
    res = apply(cbind(cell_type_real_unknown,cell_type_predict_unknown), 1, function(x) as.character(x[1]) == as.character(x[2]))
    N = length(cell_type_real_unknown)
    return(matrix(c(sum(res)/N, sum(res), N), nrow=1, ncol=3))
    
}
            
                
summary_dincta_Phi_C_accuracy <- function(meta_datas,cell_type="cell_type", 
                                      fake_cell_type = "fake_cell_type", cell_type_unknown = "unknown"){
    vals = Reduce(rbind, lapply(meta_datas, function(meta_data) 
        get_dincta_Phi_C_accuracy(meta_data,cell_type, fake_cell_type, cell_type_unknown)))
    colnames(vals) = c("accuracy","the nubmer of correct predict cell type", "total number of cells with unknow cell types")
    summary_table = data.frame(
        type=names(meta_datas)
        )
    summary_table = cbind(summary_table,vals)
    return(summary_table)
}                  
     
get_loss <- function(obj,is_dincta=TRUE){
    if(is_dincta){
#         names = c("objective_kmeans_total","objective_kmeans", "objective_kmeans_dist",
#                 "objective_kmeans_entropy",
#                "objective_kmeans_cross", "objective_kmeans_fro_cluster_error_loss",
#                "objective_kmeans_fro_R_S_loss", "objective_kmeans_Q_S_moe_loss",
#                "objective_kmeans_Q_R_loss")
        names = c("objective_kmeans", "objective_kmeans_dist",
                "objective_kmeans_entropy",
               "objective_kmeans_cross", "objective_kmeans_kl_cell_type_loss"
              )
    }else{
        names = c("objective_kmeans", "objective_kmeans_dist",
                "objective_kmeans_entropy",
               "objective_kmeans_cross")
    }
    
    loss = list()
    for(name in names){
        loss[[name]] = obj[[name]]
    }
    return(data.frame(loss))
    
}   

find_num_cell_types = function(R,eigval_threshold){
    A = diag(1/sqrt(rowSums(R))) %*% R %*% t(R) %*% diag(1/sqrt(rowSums(R))) 
    svdA = svd(A)
    d = svdA$d
    num_cell_types = sum(d>eigval_threshold)
    res= list(A=A,num_cell_types=num_cell_types)
    return(res)
}


find_gather_U <- function(A,K_sum, C=NULL, n.cell_type.min.cells =10, strict = FALSE, epsilon=1e-6,verbose = FALSE ){
    K = dim(A)[1]
    U = diag(1,K,K)
    E = A - diag(1, K,K)
    D = t(E) %*% E
    err = diag(D)
    O = D
    diag(O) = Inf
    loss = c()
    loss = append(loss,sum(err))
    i = 1
    if(verbose){
        print(paste("number of cell types: ", K, "err:", loss[1]))
    }
    
    C_find = K
    C_set = TRUE
    if(is.null(C)){
        C_set = FALSE
        C=1
        strict=FALSE
    }
    for(C_ in (K-1):C){
        join_index = arrayInd(which(O == min(O)),dim(O))
        n_max =dim(join_index)[1]
        if(n_max>1){
            err_join = rep(-Inf,n_max)
            for(ci in 1:n_max){
                irow = join_index[ci,1]
                icol = join_index[ci,2]
                if(irow<icol){
                    err_join[ci] = err[irow]+ err[icol]
                }
            }
            select_c = which.max(err_join)
            join_index = join_index[select_c,]
        }
        c1 = min(join_index)
        c2 = max(join_index)
        
        # combine D to c \times c
        index = rep(T,C_+1)
        index[c1] = F
        index[c2] = F
        index_v =which(index)
        
        D_ = matrix(Inf, C_, C_ )
        if(C_>1){
            D_[1:(C_-1),1:(C_-1)] = D[index_v,index_v]
            D_[1:(C_-1),C_] = D[index_v,c1] +  D[index_v,c2]
            D_[C_,1:(C_-1)] = D[c1,index_v] +  D[c2,index_v]
        }
        D_[C_,C_] = D[c1,c1] + D[c2,c2] + D[c1,c2] + D[c2,c1]
        
        err = diag(D_)
        err_current = sum(err)
        if(((err_current>loss[i]) | abs(err_current - loss[i])<epsilon)&(!strict)){
            if(verbose){
                 print(paste("Final, find number of cell types: ", C_find, "err:", loss[i-1]))  
                 if(C_set){
                     print(paste("If you want to reduce to the number of cell types: ",C, ", set the parameter strict=TRUE." ))
                 }
            }
           
            break
        }
        loss = append(loss,sum(err_current))
        i=i+1
        D = D_
        O = D
        diag(O) = Inf
        if(verbose){
             print(paste("number of cell types: ", C_, "err:", loss[i]))
        }
        
        u_combine =  U[c1,] + U[c2,]
        if(C_>1){
            U[1:(C_-1),] = U[index_v,]
        }
        U[C_,] =  u_combine 
        C_find = C_ 
         
#         print(paste("U:"))
#         for(ic in 1:C_ ){
#             print(paste(U[ic,], collapse = ', '))
#         }
    }
    C_sum = U[1:C_find,] %*% K_sum
    # D c_find x c_find
    # shrink the cell type until C_sum > n.cell_type.min.cells
    while((sum(C_sum < n.cell_type.min.cells)>0)&(C_find>1)){
        index = which(C_sum  < n.cell_type.min.cells)
        if(length(index)>1){
            # find the largeest error in D
            d = diag(D)
            index_ = which(d[index] == max(d[index]))
            index = index[index_] 
            index = index[1] # the cell type selected
        }
        # find which cell type to join
        o = D[index,]
        o[index] = Inf
        indexp = which(o ==min(o))
        if(length(indexp)>1){
            # find the largest error in d
            indexp_ = which(d[indexp] == max(d[indexp]))
            indexp = indexp[indexp_]
            indexp = indexp[1]
        }
        c1 = min(c(index,indexp))
        c2 = max(c(index,indexp))
        C_= C_find-1
        # join index, indxep
        # combine D to c \times c
        index = rep(T,C_+1)
        index[c1] = F
        index[c2] = F
        index_v =which(index)
        
        D_ = matrix(Inf, C_, C_ )
        if(C_>1){
            D_[1:(C_-1),1:(C_-1)] = D[index_v,index_v]
            D_[1:(C_-1),C_] = D[index_v,c1] +  D[index_v,c2]
            D_[C_,1:(C_-1)] = D[c1,index_v] +  D[c2,index_v]
        }
        D_[C_,C_] = D[c1,c1] + D[c2,c2] + D[c1,c2] + D[c2,c1]
        
        err = diag(D_)
        err_current = sum(err)
        loss = append(loss,sum(err_current))
        i=i+1
        D = D_
        O = D
        diag(O) = Inf
        if(verbose){
             print(paste("number of cell types: ", C_, "err:", loss[i]))
        }
        
        u_combine =  U[c1,] + U[c2,]
        if(C_>1){
            U[1:(C_-1),] = U[index_v,]
        }
        U[C_,] =  u_combine 
        C_find = C_ 
        C_sum = U[1:C_find,] %*% K_sum
#         print(paste("C_sum :", paste(C_sum, collpase=',')))
         
    }
    return(U[1:C_find,])
}
                                
# align the predict cell type to the known cell type 
align_cell_types <- function(Phi_C, Psi_C, verbose=F){
    C_reference = dim(Phi_C)[1]
    C_predict  = dim(Psi_C)[1]
    O = Phi_C %*% t(Psi_C)
    M = O  %*% diag(1./(colSums(O)+1e-8)) 
    # how to handle the tie
    reference_mapped = max.col(t(M))
    if(verbose){
        for(cp in 1:C_predict){
           cpm = reference_mapped[cp]
            num_max = sum(M[,cp] == M[cpm,cp])
           if(num_max > 1){
             print(paste("There is a tie for predict type ", cp, " chose to map to ", cpm))   
           }
        }
    }
    cell_types = rownames(Phi_C)
    cell_types_predict = rownames(Psi_C)
    return(list(reference_mapped, cell_types, cell_types_predict))
}
                                
infer_accuracy_core <-function(meta_data, 
                          reference_mapped,
                          cell_types,
                          cell_types_predict,
                          cell_type_predict_name = "cell_type_predict",
                          cell_type_reference_name = "cell_type"){
    cell_type_predict = as.character(meta_data[[cell_type_predict_name]])
    cell_type_reference = as.character(meta_data[[cell_type_reference_name]])
    N = length(cell_type_predict)
    cell_type_predict_mapped = rep('', N)
    N_correct = 0
    for(i in 1:N){
        ctp = cell_type_predict[i]
        ind_ctp = which(cell_types_predict == ctp)
        ind_ct= reference_mapped[ind_ctp]
        cell_type_predict_mapped[i] = cell_types[ind_ct]
        if(cell_types[ind_ct] == cell_type_reference[i]){
            N_correct = N_correct + 1
        }
    }
    return(c(N_correct/N, N_correct,N))  
} 

infer_accuracy <- function(meta_data, Psi_C, Phi_C, cell_type_predict_name = "cell_type_predict",
                          cell_type_reference_name = "cell_type", verbose=T){
    align_res =  align_cell_types(Phi_C, Psi_C, verbose)
    reference_mapped = align_res[[1]]
    cell_types = align_res[[2]]
    cell_types_predict = align_res[[3]]
    infer_res = infer_accuracy_core(meta_data, 
                          reference_mapped,
                          cell_types,
                          cell_types_predict,
                          cell_type_predict_name,
                          cell_type_reference_name)
    return(infer_res)
} 

infer_accuracy_R <- function(R,Phi_C,eigval_threshold =0.9,n.cell_type.min.cells =10,fix_C = TRUE, name=""){
    C = dim(Phi_C)[1]
    res = find_num_cell_types(R,eigval_threshold)
    K_sum = rowSums(R)
    if(fix_C){
       C_ = C 
    }else{
       C_ = res$num_cell_types
    }
    if(C!=C_){
        print(paste(name, ": C is ",C, "C_ is ", C_))
    }
        
    U = find_gather_U(res$A, K_sum,C = C_, n.cell_type.min.cells=n.cell_type.min.cells,strict = T)
    C_new = dim(U)[1] 
    if(C_new != C_){
        print(paste(name, ":C_new is ",C_new, "C_ is ", C_))
    }
    
    Psi_C = U %*% R
    new_cell_types_name = paste("new_cell_type_", 1:C_new, sep="")
    rownames(Psi_C) = new_cell_types_name
    cell_type_predict_index = apply(Psi_C, 2, which.max)
    cell_type_predict = as.character(lapply(cell_type_predict_index, function(cti){
        as.character(new_cell_types_name[cti])
    }))
    
    reference_cell_types_name = rownames(Phi_C)
    cell_type_reference_index = apply(Phi_C, 2, which.max)
    cell_type_reference = as.character(lapply(cell_type_reference_index, function(cti){
        as.character(reference_cell_types_name[cti])
    }))
    
    
    meta_data[['cell_type_predict']] 
    meta_data = data.frame(cell_type= cell_type_reference, cell_type_predict=cell_type_predict )
    infer_res = infer_accuracy(meta_data, Psi_C,Phi_C)
    res = list(infer_res=infer_res,Phi_C=Phi_C,Psi_C=Psi_C,meta_data = meta_data)
    return(res)
}   
                                
infer_accuracy_R_summary <- function(objs,Phi_C,eigval_threshold=0.9,n.cell_type.min.cells =10,fix_C = FALSE){
    names = names(objs)
    vals = Reduce(rbind,lapply(names, function(name){
        obj = objs[[name]]
        res = infer_accuracy_R(obj$R,Phi_C,eigval_threshold,n.cell_type.min.cells,fix_C,name)
        infer_res  = res$infer_res
        return(infer_res)
    }))
    colnames(vals) = c('accuracy','number of correct cells','number of total cells' ) 
    rownames(vals) = names              
    return(data.frame(vals))           
}                                