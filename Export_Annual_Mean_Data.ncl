load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

yr_len         = 150       ;
mon_len        = yr_len*12  ; 
; varlist        = (/"rsdt","rsut","rlut","rsds","rsus","rlds","rlus","hfls","hfss","prsn","rsutcs","rlutcs","rsdscs","rsuscs","rldscs"/)
varlist        = (/"clt"/)
model          = (/"ACCESS1.0","ACCESS1.3","BCC-CSM1.1","BCC-CSM1.1m","CanESM2",\
                  "CCSM4","CNRM-CM5","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
                  "GISS-E2-H","GISS-E2-R","INM-CM4","IPSL-CM5A-LR","IPSL-CM5B-LR",\
                  "MIROC5","MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-CGCM3",\
                  "NorESM1-M"/)
                              
pc             = (/"piControl"/)
ex             = (/"abrupt4xCO2"/)
m              = 0         ; 0 means 1st model
do m = 0,dimsizes(model)-1,1
    locate_pc      = ("/data/model_output/cmip5/" + model(m) + "/" + pc + "/" + varlist + "/" + "*_r1i1p1_*.nc")
    locate_ex      = ("/data/model_output/cmip5/" + model(m) + "/" + ex + "/" + varlist + "/" + "*_r1i1p1_*.nc")

    one_file_ls_pc    = systemfunc("ls " + locate_pc(0))                     ; all files of varlist(0)   
    one_file_pc       = addfile(one_file_ls_pc(0), "r")                      ; the first file of varlist(0)
    lat               = one_file_pc->lat
    lon               = one_file_pc->lon
    lat_bnds          = one_file_pc->lat_bnds
    lon_bnds          = one_file_pc->lon_bnds

    dimlat            = dimsizes(lat)
    dimlon            = dimsizes(lon)

    allvar_pc = new((/dimsizes(varlist), mon_len, dimlat, dimlon/), double, 0)
    ivar = 0
    do ivar = 0,dimsizes(varlist)-1,1
        all_file_pc                    = systemfunc("ls " + locate_pc(ivar))
        dim_file_pc                    = dimsizes(all_file_pc)
        
        ifile = 0
        dim_t_pc     = 0
        read_file_pc = 0
        var_temp_pc  = 0
        do while( dim_t_pc .lt. mon_len)
            delete([/read_file_pc, var_temp_pc/])
            read_file_pc                 = addfiles(all_file_pc( dim_file_pc-ifile-1 : dim_file_pc-1 ), "r")
            ListSetType(read_file_pc, "cat")
            var_temp_pc                  = read_file_pc[:]->$varlist(ivar)$
            dim_var_temp_pc              = dimsizes(var_temp_pc)   
            dim_t_pc                     = dim_var_temp_pc(0)                
            ifile = ifile + 1
        end do
        ; printVarSummary(read_file_pc)
        print("*** " + model(m) + " --- " + pc + " --- " + varlist(ivar) + " --- is reading " + ifile + " files... ***")
        var_pc                       = var_temp_pc(dim_t_pc-mon_len:, :, :)
        allvar_pc(ivar,:,:,:)        = var_pc    
        
        delete([/all_file_pc, dim_file_pc, read_file_pc, var_temp_pc, dim_var_temp_pc, dim_t_pc, var_pc/])
    end do
    allvar_pc!0    = "var"
    allvar_pc&var = ispan(1,dimsizes(varlist),1)
    allvar_ex = new((/dimsizes(varlist), mon_len, dimlat, dimlon/), double, 0)
    ivar = 0
    do ivar = 0,dimsizes(varlist)-1,1
        all_file_ex                    = systemfunc("ls " + locate_ex(ivar))
        dim_file_ex                    = dimsizes(all_file_ex)
        
        ifile = 0
        dim_t_ex     = 0
        read_file_ex = 0
        var_temp_ex  = 0
        do while( dim_t_ex .lt. mon_len)
            delete([/read_file_ex, var_temp_ex/])
            read_file_ex                 = addfiles(all_file_ex( 0:ifile ), "r")
            ListSetType(read_file_ex, "cat")
            var_temp_ex                  = read_file_ex[:]->$varlist(ivar)$
            dim_var_temp_ex              = dimsizes(var_temp_ex)   
            dim_t_ex                     = dim_var_temp_ex(0)                
            ifile = ifile + 1
        end do
        ; printVarSummary(read_file_ex)
        print("*** " + model(m) + " --- " + ex + " --- " + varlist(ivar) + " --- is reading " + ifile + " files... ***")
        var_ex                       = var_temp_ex(0:mon_len-1, :, :)
        allvar_ex(ivar,:,:,:)        = var_ex        
        delete([/all_file_ex, dim_file_ex, read_file_ex, var_temp_ex, dim_var_temp_ex, dim_t_ex, var_ex/])
    end do
    allvar_ex!0    = "var"
    allvar_ex&var  = ispan(1,dimsizes(varlist),1)
    delete([/locate_pc, locate_ex, one_file_ls_pc, one_file_pc, lat, lon/])
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; delete_VarAtts(allvar_pc, -1)
    ; delete_VarAtts(allvar_ex, -1)
    ; allvar_pc@_FillValue   = default_fillvalue(typeof(allvar_pc))
    ; allvar_ex@_FillValue   = default_fillvalue(typeof(allvar_ex))
    
    
    allvar_ex_an_ttemp = reshape(allvar_ex,(/dimsizes(varlist), yr_len, 12, dimlat, dimlon/))
    allvar_ex_an_temp  = dim_sum_n_Wrap(allvar_ex_an_ttemp,2) / 12.0
    copy_VarCoords(allvar_ex(:,0,:,:), allvar_ex_an_temp(:,0,:,:))
    allvar_ex_an_temp!1    = "year"
    allvar_ex_an_temp&year = ispan(1,yr_len,1)    
    
    allvar_ex_an   = allvar_ex_an_temp(year|:, var|:, lat|:, lon|:)

    allvar_cl_an   = dim_sum_n_Wrap(allvar_pc,1) / todouble(mon_len)
    copy_VarCoords(allvar_ex_an(0,:,:,:), allvar_cl_an)
    ; season = (/"DJF","JJA"/)
    ; allvar_ex_ss_temp      = month_to_seasonN(allvar_ex(time|:, var|:, lat|:, lon|:), season)
    ; allvar_ex_ss_temp!1    = "year"
    ; allvar_ex_ss_temp&year = ispan(1,yr_len,1)

    ; allvar_pc_ss           = month_to_seasonN(allvar_pc(time|:, var|:, lat|:, lon|:), season)
    ; allvar_cl_ss_temp      = dim_avg_n_Wrap(allvar_pc_ss,1)


    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; ; what the hell that i can't save it with the original "season" coordinate 
    ; allvar_ex_ss = ( allvar_ex_ss_temp + allvar_ex_ss_temp ) / 2.0
    ; copy_VarCoords(allvar_ex_an, allvar_ex_ss(0,:,:,:,:))
    ; allvar_ex_ss!0 = "season"
    ; allvar_ex_ss&season = (/0,1/)
    ; allvar_cl_ss = ( allvar_cl_ss_temp + allvar_cl_ss_temp ) / 2.0
    ; copy_VarCoords(allvar_cl_an, allvar_cl_ss(0,:,:,:))
    ; allvar_cl_ss!0 = "season"
    ; allvar_cl_ss&season = (/0,1/)
    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ivar = 0
    do ivar = 0,dimsizes(varlist)-1,1
        output     = (/varlist(ivar)+"_cl_an", varlist(ivar)+"_cl_ss", varlist(ivar)+"_ex_an", varlist(ivar)+"_ex_ss"/)
        ncFileName = "/home/sabrina/netCDF_export/singlevar/" + varlist(ivar) + "/" + varlist(ivar) + "_ex_an_" + model(m)
        system("/bin/rm -f "+ncFileName+".nc")
        print("Creating "+ncFileName+".nc ......")
        ncdf  = addfile(ncFileName + ".nc", "c")
        ncdf->lon_bnds    = lon_bnds
        ncdf->lat_bnds    = lat_bnds
        ncdf->$output(0)$ = allvar_cl_an(ivar,:,:)
        ; ncdf->$output(1)$ = allvar_cl_ss(:,ivar,:,:)
        ncdf->$output(2)$ = allvar_ex_an(:,ivar,:,:)
        ; ncdf->$output(3)$ = allvar_ex_ss(:,:,ivar,:,:)
    end do
    delete([/allvar_ex_an_ttemp, allvar_ex_an_temp, allvar_ex_an, allvar_cl_an/])
    ; delete([/allvar_ex_an_ttemp, allvar_ex_an_temp, allvar_ex_an, allvar_cl_an, allvar_ex_ss_temp, allvar_pc_ss, allvar_cl_ss_temp, allvar_ex_ss, allvar_cl_ss/])
    delete([/lon_bnds, lat_bnds, allvar_pc, allvar_ex/])
end do


end
