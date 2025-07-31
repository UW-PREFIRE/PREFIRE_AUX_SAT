from .create_aux_sat_product import create_aux_sat_product, parse_args

auxmet_fpath, ancillary_dir, input_products_dir, \
    L3_geoloc_dir, output_dir, tmpfiles_dir, atrack_idx_range = parse_args()
create_aux_sat_product(auxmet_fpath, ancillary_dir, input_products_dir, \
                       L3_geoloc_dir, output_dir, tmpfiles_dir, 
                       atrack_range_to_process=atrack_idx_range)