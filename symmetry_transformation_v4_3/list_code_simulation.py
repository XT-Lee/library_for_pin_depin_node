def pin_hex_to_hp_csv_cn3():
    import symmetry_transformation_v4_3.simulation_controller as sc
    hp = sc.simulation_controller_honeycomb_part_traps()
    prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
    output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_5973.csv'
    hp.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(output_file_csv)
    import symmetry_transformation_v4_3.list_code_analysis as lca
    sgf = lca.analyze_a_series_of_gsd_file()
    list_lcr_k_cn3 = sgf.get_cn3s_from_mysql_honeycomb_part()

def pin_hex_to_hp_fill_workflow():
    import symmetry_transformation_v4_3.simulation_controller as sc
    hpt = sc.simulation_controller_honeycomb_part_traps()
    hpt.generate_simu_index_csv_lcr081()
    prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
    index1 = 6013
    output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_'+str(int(index1))+'.csv'
    hpt.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(output_file_csv)
    
def pin_hex_to_hp_fill_brownian_workflow():
    import symmetry_transformation_v4_3.simulation_controller as sc
    hpt = sc.simulation_controller_honeycomb_part_traps()
    hpt.generate_simu_index_csv_lcr081_brownian()
    prefix_write = '/media/remote/32E2D4CCE2D49607/file_lxt/record_results_v430/honeycomb_part_pin/'
    index1 = 6013
    output_file_csv = prefix_write + 'pin_hex_to_honeycomb_part_klt_2m_gauss_b_'+str(int(index1))+'.csv'
    hpt.generate_initial_state_hexagonal_particle_honeycomb_part_trap_scan_csv(output_file_csv)
    import symmetry_transformation_v4_3.list_code_analysis as lca
    sgf = lca.analyze_a_series_of_gsd_file()
    list_lcr_k_cn3 = sgf.get_cn3s_from_csv_honeycomb_part(output_file_csv)
    print(list_lcr_k_cn3)