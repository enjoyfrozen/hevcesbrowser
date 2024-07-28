#ifndef HEVC_H_
#define HEVC_H_

#include <memory>
#include <vector>
#include <array>
#include <cstdint>
#include <cstddef>
#include <string>

namespace HEVC
{

  enum NALUnitType
  {
    NAL_TRAIL_N    = 0,
    NAL_TRAIL_R    = 1,
    NAL_TSA_N      = 2,
    NAL_TSA_R      = 3,
    NAL_STSA_N     = 4,
    NAL_STSA_R     = 5,
    NAL_RADL_N     = 6,
    NAL_RADL_R     = 7,
    NAL_RASL_N     = 8,
    NAL_RASL_R     = 9,
    NAL_BLA_W_LP   = 16,
    NAL_BLA_W_RADL = 17,
    NAL_BLA_N_LP   = 18,
    NAL_IDR_W_RADL = 19,
    NAL_IDR_N_LP   = 20,
    NAL_CRA_NUT    = 21,
    NAL_IRAP_VCL23 = 23,
    NAL_VPS        = 32,
    NAL_SPS        = 33,
    NAL_PPS        = 34,
    NAL_AUD        = 35,
    NAL_EOS_NUT    = 36,
    NAL_EOB_NUT    = 37,
    NAL_FD_NUT     = 38,
    NAL_SEI_PREFIX = 39,
    NAL_SEI_SUFFIX = 40,
  };

  class NALHeader
  {
  public:
    NALUnitType           type;
    uint8_t               layer_id;
    uint8_t               temporal_id_plus1;

    bool operator == (const NALHeader &) const;
  };

  class ProfileTierLevel
  {
  public:
    uint8_t                general_profile_space;
    uint8_t                general_tier_flag;
    uint8_t                general_profile_idc;
    uint8_t                general_profile_compatibility_flag[32];
    uint8_t                general_progressive_source_flag;
    uint8_t                general_interlaced_source_flag;
    uint8_t                general_non_packed_constraint_flag;
    uint8_t                general_frame_only_constraint_flag;

    uint8_t general_max_12bit_constraint_flag;
    uint8_t general_max_10bit_constraint_flag;
    uint8_t general_max_8bit_constraint_flag;
    uint8_t general_max_422chroma_constraint_flag;
    uint8_t general_max_420chroma_constraint_flag;
    uint8_t general_max_monochrome_constraint_flag;
    uint8_t general_intra_constraint_flag;
    uint8_t general_one_picture_only_constraint_flag;
    uint8_t general_lower_bit_rate_constraint_flag;
    uint8_t general_inbld_flag;

    uint8_t                general_level_idc;
    std::vector<uint8_t>   sub_layer_profile_present_flag;
    std::vector<uint8_t>   sub_layer_level_present_flag;
    std::vector<uint8_t>   sub_layer_profile_space;
    std::vector<uint8_t>   sub_layer_tier_flag;
    std::vector<uint8_t>   sub_layer_profile_idc;
    std::vector< std::vector< uint8_t> >
                           sub_layer_profile_compatibility_flag;
    std::vector<uint8_t>   sub_layer_progressive_source_flag;
    std::vector<uint8_t>   sub_layer_interlaced_source_flag;
    std::vector<uint8_t>   sub_layer_non_packed_constraint_flag;
    std::vector<uint8_t>   sub_layer_frame_only_constraint_flag;
    std::vector<uint8_t>   sub_layer_level_idc;


    std::vector<uint8_t>  sub_layer_max_12bit_constraint_flag;
    std::vector<uint8_t>  sub_layer_max_10bit_constraint_flag;
    std::vector<uint8_t>  sub_layer_max_8bit_constraint_flag;
    std::vector<uint8_t>  sub_layer_max_422chroma_constraint_flag;
    std::vector<uint8_t>  sub_layer_max_420chroma_constraint_flag;
    std::vector<uint8_t>  sub_layer_max_monochrome_constraint_flag;
    std::vector<uint8_t>  sub_layer_intra_constraint_flag;
    std::vector<uint8_t>  sub_layer_one_picture_only_constraint_flag;
    std::vector<uint8_t>  sub_layer_lower_bit_rate_constraint_flag;
    std::vector<uint8_t>  sub_layer_inbld_flag;

    void toDefault();

    bool operator == (const ProfileTierLevel &) const;
  };

  class SubLayerHrdParameters
  {
  public:
    std::vector<uint32_t>       bit_rate_value_minus1;
    std::vector<uint32_t>       cpb_size_value_minus1;
    std::vector<uint32_t>       cpb_size_du_value_minus1;
    std::vector<uint32_t>       bit_rate_du_value_minus1;
    std::vector<uint8_t>        cbr_flag;

    void toDefault();

    bool operator == (const SubLayerHrdParameters &) const;
  };


  class ScalingListData
  {
  public:
    std::vector< std::vector< uint8_t> >    scaling_list_pred_mode_flag;
    std::vector< std::vector< uint32_t> >   scaling_list_pred_matrix_id_delta;
    std::vector< std::vector< uint32_t> >   scaling_list_dc_coef_minus8;
    std::vector<std::vector< std::vector< uint32_t> > >
                                            scaling_list_delta_coef;

    void toDefault();

    bool operator == (const ScalingListData &) const;
  };

  class HrdParameters
  {
  public:
    uint8_t               nal_hrd_parameters_present_flag;
    uint8_t               vcl_hrd_parameters_present_flag;
    uint8_t               sub_pic_hrd_params_present_flag;
    uint8_t               tick_divisor_minus2;
    uint8_t               du_cpb_removal_delay_increment_length_minus1;
    uint8_t               sub_pic_cpb_params_in_pic_timing_sei_flag;
    uint8_t               dpb_output_delay_du_length_minus1;
    uint8_t               bit_rate_scale;
    uint8_t               cpb_size_scale;
    uint8_t               cpb_size_du_scale;
    uint8_t               initial_cpb_removal_delay_length_minus1;
    uint8_t               au_cpb_removal_delay_length_minus1;
    uint8_t               dpb_output_delay_length_minus1;
    std::vector<uint8_t>  fixed_pic_rate_general_flag;
    std::vector<uint8_t>  fixed_pic_rate_within_cvs_flag;
    std::vector<uint32_t> elemental_duration_in_tc_minus1;
    std::vector<uint8_t>  low_delay_hrd_flag;
    std::vector<uint32_t> cpb_cnt_minus1;
    std::vector<SubLayerHrdParameters>
                          nal_sub_layer_hrd_parameters;
    std::vector<SubLayerHrdParameters>
                          vcl_sub_layer_hrd_parameters;

    void toDefault();

    bool operator == (const HrdParameters &) const;
 };

  class VideoSignalInfo
  {
  public:
      uint8_t video_vps_format;
      uint8_t video_full_range_vps_flag;
      uint8_t colour_primaries_vps;
      uint8_t transfer_characteristics_vps;
      uint8_t matrix_coeffs_vps;
      
      void toDefault();

      bool operator == (const VideoSignalInfo&) const;
  };


  //vps_vui_bsp_hrd_params
  class VpsVuiBspHrdParams {
  public:
      uint32_t vps_num_add_hrd_params;
      std::vector<uint8_t> cprms_add_present_flag;
      std::vector<uint32_t> num_sub_layer_hrd_minus1;
      std::vector < HrdParameters> hrd_parameters;
      std::vector<uint32_t> num_signalled_partitioning_schemes;
      std::vector < std::vector<uint32_t>> num_partitions_in_scheme_minus1;
      std::vector < std::vector < std::vector < std::vector < uint8_t> > >> layer_included_in_partition_flag;
      std::vector < std::vector < std::vector < uint32_t> > > num_bsp_schedules_minus1;
      std::vector < std::vector < std::vector < std::vector < std::vector < uint32_t> > >>> bsp_hrd_idx;
      std::vector < std::vector < std::vector < std::vector < std::vector < uint32_t> > >>> bsp_sched_idx;

      void toDefault();

  };


  //zorro add VpsVui
  class VpsVui
  {
  public:
      uint8_t cross_layer_pic_type_aligned_flag;
      uint8_t cross_layer_irap_aligned_flag;
      uint8_t all_layers_idr_aligned_flag;
      uint8_t bit_rate_present_vps_flag;
      uint8_t pic_rate_present_vps_flag;
      std::vector< std::vector< uint8_t> > bit_rate_present_flag;
      std::vector< std::vector< uint8_t> >  pic_rate_present_flag;
      std::vector< std::vector< uint16_t> > avg_bit_rate;
      std::vector< std::vector< uint16_t> > max_bit_rate;
      std::vector< std::vector< uint8_t> > constant_pic_rate_idc;
      std::vector< std::vector< uint16_t> > avg_pic_rate;

      uint8_t video_signal_info_idx_present_flag;
      uint8_t vps_num_video_signal_info_minus1;

      std::vector < VideoSignalInfo> video_signal_info;
      std::vector< uint8_t> vps_video_signal_info_idx;
      uint8_t tiles_not_in_use_flag;
      std::vector< uint8_t> tiles_in_use_flag;
      std::vector< uint8_t> loop_filter_not_across_tiles_flag;
      std::vector < std::vector< uint8_t>> tile_boundaries_aligned_flag;

      uint8_t wpp_not_in_use_flag;
      std::vector< uint8_t>  wpp_in_use_flag;
      uint8_t single_layer_for_non_irap_flag;
      uint8_t higher_layer_irap_skip_flag;
      uint8_t ilp_restricted_ref_layers_flag;
      std::vector < std::vector< uint8_t>> min_spatial_segment_offset_plus1;
      std::vector < std::vector< uint8_t>> ctu_based_offset_enabled_flag;
      std::vector < std::vector< uint8_t>> min_horizontal_ctu_offset_plus1;
      uint8_t vps_vui_bsp_hrd_present_flag;

      VpsVuiBspHrdParams vps_vui_bsp_hrd_params;
      std::vector< uint8_t> base_layer_parameter_set_compatibility_flag;

      void toDefault();

  };

  //dpb_size
  class DpbSize
  {
  public:
      std::vector< uint8_t> sub_layer_flag_info_present_flag;
      std::vector < std::vector< uint8_t>> sub_layer_dpb_info_present_flag;
      std::vector<std::vector< std::vector< uint8_t> > > max_vps_dec_pic_buffering_minus1;
      std::vector< std::vector< uint8_t> > max_vps_num_reorder_pics;
      std::vector< std::vector< uint8_t> > max_vps_latency_increase_plus1;

      void toDefault();

      bool operator == (const DpbSize&) const;
  };

  //rep_format
  class RepFormat
  {
  public:
      uint16_t pic_width_vps_in_luma_samples;
      uint16_t pic_height_vps_in_luma_samples;
      uint8_t chroma_and_bit_depth_vps_present_flag;
      uint8_t chroma_format_vps_idc;
      uint8_t separate_colour_plane_vps_flag;
      uint8_t bit_depth_vps_luma_minus8;
      uint8_t bit_depth_vps_chroma_minus8;
      uint8_t conformance_window_vps_flag;
      uint32_t conf_win_vps_left_offset;
      uint32_t conf_win_vps_right_offset;
      uint32_t conf_win_vps_top_offset;
      uint32_t conf_win_vps_bottom_offset;


      void toDefault();

      bool operator == (const RepFormat&) const;
  };

  //zorro add  vps_extension
  class VpsExtension
  {
  public:
    //profile_tier_level( 0, vps_max_sub_layers_minus1 )

    uint8_t splitting_flag;
    std::vector<uint8_t> scalability_mask_flag;

    std::vector<uint8_t>  dimension_id_len_minus1;
    uint8_t vps_nuh_layer_id_present_flag;
    std::vector<uint8_t> layer_id_in_nuh;
    std::vector< std::vector< uint8_t> > dimension_id;
    uint8_t view_id_len;
    std::vector<uint8_t>  view_id_val;
    std::vector< std::vector< uint8_t> > direct_dependency_flag;
    uint32_t num_add_layer_sets;
    std::vector< std::vector< uint8_t> > highest_layer_idx_plus1;
    uint8_t vps_sub_layers_max_minus1_present_flag;
    std::vector< uint8_t> sub_layers_vps_max_minus1;
    uint8_t max_tid_ref_present_flag;
    std::vector< std::vector< uint8_t> > max_tid_il_ref_pics_plus1;
    uint8_t default_ref_layers_active_flag;
    uint8_t vps_num_profile_tier_level_minus1;
    std::vector<uint8_t>  vps_profile_present_flag;
    //profile_tier_level( vps_profile_present_flag[ i ], vps_max_sub_layers_minus1 )
    uint8_t num_add_olss;
    uint8_t default_output_layer_idc;
    //uint32_t NumOutputLayerSets() {  return num_add_olss + NumLayerSets();  }
    std::vector<uint8_t> layer_set_idx_for_ols_minus1;
    std::vector< std::vector< uint8_t> > output_layer_flag;
    std::vector< std::vector< uint8_t> > profile_tier_level_idx;
    std::vector< uint8_t> alt_output_layer_flag;
    uint8_t vps_num_rep_formats_minus1;
    std::vector < RepFormat> rep_format;
    uint8_t rep_format_idx_present_flag;
    std::vector<uint8_t> vps_rep_format_idx;
    uint8_t max_one_active_ref_layer_flag;
    uint8_t vps_poc_lsb_aligned_flag;
    std::vector<uint8_t> poc_lsb_not_present_flag;
    DpbSize dpb_size;
    uint8_t direct_dep_type_len_minus2;
    uint8_t direct_dependency_all_layers_flag;
    uint32_t direct_dependency_all_layers_type;
    std::vector< std::vector< uint8_t> >  direct_dependency_type;
    uint32_t vps_non_vui_extension_length;
    std::vector<uint8_t> vps_non_vui_extension_data_byte;
    uint8_t vps_vui_present_flag;
    //alignbyte
    VpsVui vps_vui;

    void toDefault();

    bool operator == (const VpsExtension&) const;
  };
  
  //zorro add  vps_3d_extension
  class Vps3dExtension
  {
  public:
      uint32_t cp_precision;
      std::vector<uint8_t> num_cp;
      std::vector<uint8_t> cp_in_slice_segment_header_flag;
      std::vector<std::vector< std::vector< uint32_t> > > cp_ref_voi;
      std::vector<std::vector< std::vector< uint32_t> > > vps_cp_scale;
      std::vector<std::vector< std::vector< uint32_t> > > vps_cp_off;
      std::vector<std::vector< std::vector< uint32_t> > > vps_cp_inv_scale_plus_scale;
      std::vector<std::vector< std::vector< uint32_t> > > vps_cp_inv_off_plus_off;

      void toDefault();

      bool operator == (const Vps3dExtension&) const;
  };

  //zorro add  sps_3d_extension
  class Sps3dExtension
  {
  public:
      uint8_t iv_di_mc_enabled_flag[2];
      uint8_t iv_mv_scal_enabled_flag[2];
      uint32_t log2_ivmc_sub_pb_size_minus3 [2];
      uint8_t iv_res_pred_enabled_flag[2];
      uint8_t depth_ref_enabled_flag[2];
      uint8_t vsp_mc_enabled_flag[2];
      uint8_t dbbp_enabled_flag[2];
      uint8_t tex_mc_enabled_flag[2];
      uint8_t log2_texmc_sub_pb_size_minus3[2];
      uint8_t intra_contour_enabled_flag [2];
      uint8_t intra_dc_only_wedge_enabled_flag[2];
      uint8_t cqt_cu_part_pred_enabled_flag[2];
      uint8_t inter_dc_only_enabled_flag[2];
      uint8_t skip_intra_enabled_flag [2];
      
      void toDefault();

      bool operator == (const Sps3dExtension&) const;
  };

  //sps_range_extension
  class SpsRangeExtension
  {
  public:
      uint8_t transform_skip_rotation_enabled_flag;
      uint8_t transform_skip_context_enabled_flag;
      uint8_t implicit_rdpcm_enabled_flag;
      uint8_t explicit_rdpcm_enabled_flag;
      
      uint8_t extended_precision_processing_flag;

      uint8_t intra_smoothing_disabled_flag;
      uint8_t high_precision_offsets_enabled_flag;
      uint8_t persistent_rice_adaptation_enabled_flag;
      uint8_t cabac_bypass_alignment_enabled_flag;

      void toDefault();
  };

  //pps_range_extension
  class PpsRangeExtension {
  public:
      uint32_t log2_max_transform_skip_block_size_minus2;
      uint8_t cross_component_prediction_enabled_flag;
      uint8_t chroma_qp_offset_list_enabled_flag;
      uint32_t diff_cu_chroma_qp_offset_depth;
      uint32_t chroma_qp_offset_list_len_minus1;
      std::vector<int32_t> cb_qp_offset_list;
      std::vector<int32_t> cr_qp_offset_list;
      
      uint32_t log2_sao_offset_scale_luma;
      uint32_t log2_sao_offset_scale_chroma;

      void toDefault();
  };

  //colour_mapping_table
  class ColourMappingTable
  {
  public:
      uint8_t num_cm_ref_layers_minus1;
      std::vector<uint8_t> cm_ref_layer_id;
      uint8_t cm_octant_depth;
      uint8_t cm_y_part_num_log2;
      uint8_t luma_bit_depth_cm_input_minus8;
      uint8_t chroma_bit_depth_cm_input_minus8;
      uint8_t luma_bit_depth_cm_output_minus8;
      uint8_t chroma_bit_depth_cm_output_minus8;
      uint8_t cm_res_quant_bits;
      uint8_t cm_delta_flc_bits_minus1;
      uint8_t cm_adapt_threshold_u_delta;
      uint8_t cm_adapt_threshold_v_delta;
      //colour_mapping_octants

      void toDefault();
  };

  //pps_multilayer_extension
  class PpsMultilayerExtension
  {
  public:
      uint8_t poc_reset_info_present_flag;
      uint8_t pps_infer_scaling_list_flag;
      uint8_t pps_scaling_list_ref_layer_id;
      uint16_t num_ref_loc_offsets;
      //for (i = 0; i < num_ref_loc_offsets; i++) 
      //{
          std::vector<uint8_t> ref_loc_offset_layer_id;
          std::vector<uint8_t> scaled_ref_layer_offset_present_flag;
          std::vector<uint8_t> scaled_ref_layer_left_offset;
          std::vector<uint8_t> scaled_ref_layer_top_offset;
          std::vector<uint8_t> scaled_ref_layer_right_offset;
          std::vector<uint8_t> scaled_ref_layer_bottom_offset;

      //}
      std::vector<uint8_t> ref_region_offset_present_flag;
      std::vector<uint8_t> ref_region_left_offset;
      std::vector<uint8_t> ref_region_top_offset;
      std::vector<uint8_t> ref_region_right_offset;
      std::vector<uint8_t> ref_region_bottom_offset;
      std::vector<uint8_t> resample_phase_set_present_flag;
      std::vector<uint8_t> phase_hor_luma;
      std::vector<uint8_t> phase_ver_luma;
      std::vector<uint8_t> phase_hor_chroma_plus8;
      std::vector<uint8_t> phase_ver_chroma_plus8;
      uint8_t colour_mapping_enabled_flag;
      ColourMappingTable colour_mapping_table;

      void toDefault();

  };


  //delta_dlt
  class DeltaDlt
  {
  public:
      uint32_t num_val_delta_dlt;
      uint32_t max_diff;
      uint32_t min_diff_minus1;
      uint32_t delta_dlt_val0;
      std::vector< uint8_t> delta_val_diff_minus_min;

      void toDefault();

      bool operator == (const DeltaDlt&) const;
  };

  //pps_3d_extension(
  class Pps3dExtension 
  {
  public:
      uint8_t dlts_present_flag;
      uint8_t pps_depth_layers_minus1;
      uint8_t pps_bit_depth_for_depth_layers_minus8;
      std::vector<uint8_t> dlt_flag;
      std::vector<uint8_t> dlt_pred_flag;
      std::vector<uint8_t> dlt_val_flags_present_flag;
      std::vector< std::vector< uint8_t> >  dlt_value_flag;

      std::vector < DeltaDlt> delta_dlt;

      void toDefault();

      bool operator == (const Pps3dExtension&) const;
  };

  class ShortTermRefPicSet
  {
  public:
    uint8_t                   inter_ref_pic_set_prediction_flag;
    uint32_t                  delta_idx_minus1;
    uint8_t                   delta_rps_sign;
    uint32_t                  abs_delta_rps_minus1;
    std::vector<uint8_t>      used_by_curr_pic_flag;
    std::vector<uint8_t>      use_delta_flag;
    uint32_t                  num_negative_pics;
    uint32_t                  num_positive_pics;
    std::vector<uint32_t>     delta_poc_s0_minus1;
    std::vector<uint8_t>      used_by_curr_pic_s0_flag;
    std::vector<uint32_t>     delta_poc_s1_minus1;
    std::vector<uint8_t>      used_by_curr_pic_s1_flag;

    void toDefault();

    bool operator == (const ShortTermRefPicSet &) const;
  };

  class RefPicListModification
  {
  public:
    uint8_t                ref_pic_list_modification_flag_l0;
    std::vector<uint32_t>  list_entry_l0;
    uint8_t                ref_pic_list_modification_flag_l1;
    std::vector<uint32_t>  list_entry_l1;

    void toDefault();

    bool operator == (const RefPicListModification &) const;
  };

  class VuiParameters
  {
  public:
    uint8_t          aspect_ratio_info_present_flag;
    uint8_t          aspect_ratio_idc;
    uint16_t         sar_width;
    uint16_t         sar_height;
    uint8_t          overscan_info_present_flag;
    uint8_t          overscan_appropriate_flag;
    uint8_t          video_signal_type_present_flag;
    uint8_t          video_format;
    uint8_t          video_full_range_flag;
    uint8_t          colour_description_present_flag;
    uint8_t          colour_primaries;
    uint8_t          transfer_characteristics;
    uint8_t          matrix_coeffs;
    uint8_t          chroma_loc_info_present_flag;
    uint32_t         chroma_sample_loc_type_top_field;
    uint32_t         chroma_sample_loc_type_bottom_field;
    uint8_t          neutral_chroma_indication_flag;
    uint8_t          field_seq_flag;
    uint8_t          frame_field_info_present_flag;
    uint8_t          default_display_window_flag;
    uint32_t         def_disp_win_left_offset;
    uint32_t         def_disp_win_right_offset;
    uint32_t         def_disp_win_top_offset;
    uint32_t         def_disp_win_bottom_offset;
    uint8_t          vui_timing_info_present_flag;
    uint32_t         vui_num_units_in_tick;
    uint32_t         vui_time_scale;
    uint8_t          vui_poc_proportional_to_timing_flag;
    uint32_t         vui_num_ticks_poc_diff_one_minus1;
    uint8_t          vui_hrd_parameters_present_flag;
    HrdParameters    hrd_parameters;
    uint8_t          bitstream_restriction_flag;
    uint8_t          tiles_fixed_structure_flag;
    uint8_t          motion_vectors_over_pic_boundaries_flag;
    uint8_t          restricted_ref_pic_lists_flag;
    uint32_t         min_spatial_segmentation_idc;
    uint32_t         max_bytes_per_pic_denom;
    uint32_t         max_bits_per_min_cu_denom;
    uint32_t         log2_max_mv_length_horizontal;
    uint32_t         log2_max_mv_length_vertical;

    void toDefault();

    bool operator == (const VuiParameters &) const;

  };


  class SeiPayload
  {
  public:
    virtual ~SeiPayload() {};
  };

  class SeiMessage
  {
  public:
    enum PayloadType
    {
        BUFFERING_PERIOD                     = 0,
        PICTURE_TIMING                       = 1,
        PAN_SCAN_RECT                        = 2,
        FILLER_PAYLOAD                       = 3,
        USER_DATA_REGISTERED_ITU_T_T35       = 4,
        USER_DATA_UNREGISTERED               = 5,
        RECOVERY_POINT                       = 6,
        SCENE_INFO                           = 9,
        FULL_FRAME_SNAPSHOT                  = 15,
        PROGRESSIVE_REFINEMENT_SEGMENT_START = 16,
        PROGRESSIVE_REFINEMENT_SEGMENT_END   = 17,
        FILM_GRAIN_CHARACTERISTICS           = 19,
        POST_FILTER_HINT                     = 22,
        TONE_MAPPING_INFO                    = 23,
        FRAME_PACKING                        = 45,
        DISPLAY_ORIENTATION                  = 47,
        SOP_DESCRIPTION                      = 128,
        ACTIVE_PARAMETER_SETS                = 129,
        DECODING_UNIT_INFO                   = 130,
        TEMPORAL_LEVEL0_INDEX                = 131,
        DECODED_PICTURE_HASH                 = 132,
        SCALABLE_NESTING                     = 133,
        REGION_REFRESH_INFO                  = 134,
        NO_DISPLAY                           = 135,
        TIME_CODE                            = 136,
        MASTERING_DISPLAY_INFO               = 137,
        SEGM_RECT_FRAME_PACKING              = 138,
        TEMP_MOTION_CONSTRAINED_TILE_SETS    = 139,
        CHROMA_RESAMPLING_FILTER_HINT        = 140,
        KNEE_FUNCTION_INFO                   = 141,
        COLOUR_REMAPPING_INFO                = 142,
        CONTENT_LIGHT_LEVEL_INFO             = 144,
        ALTERNATIVE_TRANSFER_CHARACTERISTICS = 147,
    };

    uint32_t           num_payload_type_ff_bytes;
    uint32_t           num_payload_size_ff_bytes;
    uint8_t            last_payload_type_byte;
    uint8_t            last_payload_size_byte;
    std::shared_ptr<SeiPayload>
                       sei_payload;

    virtual void toDefault();
  };


  class DecodedPictureHash: public SeiPayload
  {
  public:
    uint8_t                                   hash_type;
    std::vector<std::array<uint8_t, 16> >     picture_md5;
    std::vector<uint16_t>                     picture_crc;
    std::vector<uint32_t>                     picture_checksum;

    void toDefault();
  };

  class UserDataUnregistered: public SeiPayload
  {
  public:
    uint8_t                               uuid_iso_iec_11578[16];
    std::vector<uint8_t>                  user_data_payload_byte;

    void toDefault() {};
  };

  class SceneInfo: public SeiPayload
  {
  public:
    uint8_t                               scene_info_present_flag;
    uint8_t                               prev_scene_id_valid_flag;
    uint32_t                              scene_id;
    uint32_t                              scene_transition_type;
    uint32_t                              second_scene_id;

    void toDefault() {};
  };

  class FullFrameSnapshot: public SeiPayload
  {
  public:
    uint32_t                              snapshot_id;

    void toDefault() {};
  };

  class ProgressiveRefinementSegmentStart: public SeiPayload
  {
  public:
    uint32_t                              progressive_refinement_id;
    uint32_t                              pic_order_cnt_delta;

    void toDefault() {};
  };

  class ProgressiveRefinementSegmentEnd: public SeiPayload
  {
  public:
    uint32_t                              progressive_refinement_id;

    void toDefault() {};
  };

  class BufferingPeriod: public SeiPayload
  {
  public:
    uint32_t                           bp_seq_parameter_set_id;
    uint8_t                            irap_cpb_params_present_flag;
    uint32_t                           cpb_delay_offset;
    uint32_t                           dpb_delay_offset;
    uint8_t                            concatenation_flag;
    uint32_t                           au_cpb_removal_delay_delta_minus1;
    std::vector<uint32_t>              nal_initial_cpb_removal_delay;
    std::vector<uint32_t>              nal_initial_cpb_removal_offset;
    std::vector<uint32_t>              nal_initial_alt_cpb_removal_delay;
    std::vector<uint32_t>              nal_initial_alt_cpb_removal_offset;
    std::vector<uint32_t>              vcl_initial_cpb_removal_delay;
    std::vector<uint32_t>              vcl_initial_cpb_removal_offset;
    std::vector<uint32_t>              vcl_initial_alt_cpb_removal_delay;
    std::vector<uint32_t>              vcl_initial_alt_cpb_removal_offset;

    void toDefault();
  };

  class FillerPayload: public SeiPayload
  {
  public:
    void toDefault() {};
  };

  class PicTiming: public SeiPayload
  {
  public:
    uint8_t                            pic_struct;
    uint8_t                            source_scan_type;
    uint8_t                            duplicate_flag;
    uint32_t                           au_cpb_removal_delay_minus1;
    uint32_t                           pic_dpb_output_delay;
    uint32_t                           pic_dpb_output_du_delay;
    uint32_t                           num_decoding_units_minus1;
    uint8_t                            du_common_cpb_removal_delay_flag;
    uint32_t                           du_common_cpb_removal_delay_increment_minus1;
    std::vector<uint32_t>              num_nalus_in_du_minus1;
    std::vector<uint32_t>              du_cpb_removal_delay_increment_minus1;

    void toDefault();
  };

  class RecoveryPoint: public SeiPayload
  {
  public:
    uint32_t                           recovery_poc_cnt;
    uint8_t                            exact_match_flag;
    uint8_t                            broken_link_flag;
    void toDefault() {};
  };


  class ActiveParameterSets: public SeiPayload
  {
  public:
    uint8_t                 active_video_parameter_set_id;
    uint8_t                 self_contained_cvs_flag;
    uint8_t                 no_parameter_set_update_flag;
    uint32_t                num_sps_ids_minus1;
    std::vector<uint32_t>   active_seq_parameter_set_id;

    void toDefault() {};
  };


  class TemporalLevel0Index: public SeiPayload
  {
  public:
    uint8_t                 temporal_sub_layer_zero_idx;
    uint8_t                 irap_pic_id;

    void toDefault() {};
  };

  class RegionRefreshInfo: public SeiPayload
  {
  public:
    uint8_t                 refreshed_region_flag;

    void toDefault() {};
  };

  class ToneMapping: public SeiPayload
  {
  public:
    uint32_t                tone_map_id;
    uint8_t                 tone_map_cancel_flag;
    uint8_t                 tone_map_persistence_flag;
    uint8_t                 coded_data_bit_depth;
    uint8_t                 target_bit_depth;
    uint32_t                tone_map_model_id;
    uint32_t                min_value;
    uint32_t                max_value;
    uint32_t                sigmoid_midpoint;
    uint32_t                sigmoid_width;
    std::vector<uint32_t>   start_of_coded_interval;
    uint16_t                num_pivots;
    std::vector<uint32_t>   coded_pivot_value;
    std::vector<uint32_t>   target_pivot_value;
    uint8_t                 camera_iso_speed_idc;
    uint32_t                camera_iso_speed_value;
    uint8_t                 exposure_index_idc;
    uint32_t                exposure_index_value;
    uint8_t                 exposure_compensation_value_sign_flag;
    uint16_t                exposure_compensation_value_numerator;
    uint16_t                exposure_compensation_value_denom_idc;
    uint32_t                ref_screen_luminance_white;
    uint32_t                extended_range_white_level;
    uint16_t                nominal_black_level_code_value;
    uint16_t                nominal_white_level_code_value;
    uint16_t                extended_white_level_code_value;

    void toDefault() {};
  };


  class FramePacking: public SeiPayload
  {
  public:
    uint32_t                frame_packing_arrangement_id;
    uint8_t                 frame_packing_arrangement_cancel_flag;
    uint8_t                 frame_packing_arrangement_type;
    uint8_t                 quincunx_sampling_flag;
    uint8_t                 content_interpretation_type;
    uint8_t                 spatial_flipping_flag;
    uint8_t                 frame0_flipped_flag;
    uint8_t                 field_views_flag;
    uint8_t                 current_frame_is_frame0_flag;
    uint8_t                 frame0_self_contained_flag;
    uint8_t                 frame1_self_contained_flag;
    uint8_t                 frame0_grid_position_x;
    uint8_t                 frame0_grid_position_y;
    uint8_t                 frame1_grid_position_x;
    uint8_t                 frame1_grid_position_y;
    uint8_t                 frame_packing_arrangement_reserved_byte;
    uint8_t                 frame_packing_arrangement_persistence_flag;
    uint8_t                 upsampled_aspect_ratio_flag;

    void toDefault() {};
  };


  class DisplayOrientation: public SeiPayload
  {
  public:
    uint8_t                 display_orientation_cancel_flag;
    uint8_t                 hor_flip;
    uint8_t                 ver_flip;
    uint16_t                anticlockwise_rotation;
    uint8_t                 display_orientation_persistence_flag;

    void toDefault() {};
  };

  class SOPDescription: public SeiPayload
  {
  public:
    uint32_t                sop_seq_parameter_set_id;
    uint32_t                num_entries_in_sop_minus1;
    std::vector<uint8_t>    sop_vcl_nut;
    std::vector<uint8_t>    sop_temporal_id;
    std::vector<uint32_t>   sop_short_term_rps_idx;
    std::vector<uint32_t>   sop_poc_delta;

    void toDefault() {};
  };

  class TimeCode: public SeiPayload
  {
  public:
    uint8_t                       num_clock_ts;
    std::vector<uint8_t>          clock_time_stamp_flag;
    std::vector<uint8_t>          nuit_field_based_flag;
    std::vector<uint8_t>          counting_type;
    std::vector<uint8_t>          full_timestamp_flag;
    std::vector<uint8_t>          discontinuity_flag;
    std::vector<uint8_t>          cnt_dropped_flag;
    std::vector<uint16_t>         n_frames;
    std::vector<uint8_t>          seconds_value;
    std::vector<uint8_t>          minutes_value;
    std::vector<uint8_t>          hours_value;
    std::vector<uint8_t>          seconds_flag;
    std::vector<uint8_t>          minutes_flag;
    std::vector<uint8_t>          hours_flag;
    std::vector<uint8_t>          time_offset_length;
    std::vector<uint32_t>         time_offset_value;

    void toDefault() {};
  };

  class MasteringDisplayInfo: public SeiPayload
  {
  public:
    uint16_t      display_primary_x[3];
    uint16_t      display_primary_y[3];
    uint16_t      white_point_x;
    uint16_t      white_point_y;
    uint32_t      max_display_mastering_luminance;
    uint32_t      min_display_mastering_luminance;

    void toDefault();
  };

  class SegmRectFramePacking: public SeiPayload
  {
  public:
    uint8_t      segmented_rect_frame_packing_arrangement_cancel_flag;
    uint8_t      segmented_rect_content_interpretation_type;
    uint8_t      segmented_rect_frame_packing_arrangement_persistence;

    void toDefault();
  };

  class KneeFunctionInfo: public SeiPayload
  {
  public:
    uint32_t                  knee_function_id;
    uint8_t                   knee_function_cancel_flag;
    uint8_t                   knee_function_persistence_flag;
    uint32_t                  input_d_range;
    uint32_t                  input_disp_luminance;
    uint32_t                  output_d_range;
    uint32_t                  output_disp_luminance;
    uint32_t                  num_knee_points_minus1;
    std::vector<uint16_t>     input_knee_point;
    std::vector<uint16_t>     output_knee_point;

    void toDefault();
  };

  class ChromaResamplingFilterHint: public SeiPayload
  {
  public:
    uint8_t                   ver_chroma_filter_idc;
    uint8_t                   hor_chroma_filter_idc;
    uint8_t                   ver_filtering_field_processing_flag;
    uint32_t                  target_format_idc;
    uint32_t                  num_vertical_filters;
    std::vector<uint32_t>     ver_tap_length_minus_1;
    std::vector<
     std::vector<int32_t> >   ver_filter_coeff;

    uint32_t                  num_horizontal_filters;
    std::vector<uint32_t>     hor_tap_length_minus_1;
    std::vector<
     std::vector<int32_t> >   hor_filter_coeff;

    void toDefault();
  };

  class ColourRemappingInfo: public SeiPayload
  {
  public:
    uint32_t                     colour_remap_id;
    uint8_t                      colour_remap_cancel_flag;
    uint8_t                      colour_remap_persistence_flag;
    uint8_t                      colour_remap_video_signal_info_present_flag;
    uint8_t                      colour_remap_full_range_flag;
    uint8_t                      colour_remap_primaries;
    uint8_t                      colour_remap_transfer_function;
    uint8_t                      colour_remap_matrix_coefficients;
    uint8_t                      colour_remap_input_bit_depth;
    uint8_t                      colour_remap_bit_depth;
    uint8_t                      pre_lut_num_val_minus1[3];
    std::vector<uint32_t>        pre_lut_coded_value[3];
    std::vector<uint32_t>        pre_lut_target_value[3];
    uint8_t                      colour_remap_matrix_present_flag;
    uint8_t                      log2_matrix_denom;
    int32_t                      colour_remap_coeffs[3][3];
    uint8_t                      post_lut_num_val_minus1[3];
    std::vector<uint32_t>        post_lut_coded_value[3];
    std::vector<uint32_t>        post_lut_target_value[3];

    void toDefault() {};
  };

  class ContentLightLevelInfo: public SeiPayload
  {
  public:
    uint16_t      max_content_light_level;
    uint16_t      max_pic_average_light_level;

    void toDefault() {};
  };

  class AlternativeTransferCharacteristics: public SeiPayload
  {
  public:
    uint16_t      alternative_transfer_characteristics;

    void toDefault() {};
  };

  class PredWeightTable
  {
  public:
    uint32_t                  luma_log2_weight_denom;
    int32_t                   delta_chroma_log2_weight_denom;
    std::vector<uint8_t>      luma_weight_l0_flag;
    std::vector<uint8_t>      chroma_weight_l0_flag;
    std::vector<int32_t>      delta_luma_weight_l0;
    std::vector<int32_t>      luma_offset_l0;
    std::vector<std::array<int32_t, 2> >
                              delta_chroma_weight_l0;
    std::vector<std::array<int32_t, 2> >
                              delta_chroma_offset_l0;
    std::vector<uint8_t>      luma_weight_l1_flag;
    std::vector<uint8_t>      chroma_weight_l1_flag;
    std::vector<int32_t>      delta_luma_weight_l1;
    std::vector<int32_t>      luma_offset_l1;
    std::vector<std::array<int32_t, 2> >
                              delta_chroma_weight_l1;
    std::vector<std::array<int32_t, 2> >
                              delta_chroma_offset_l1;

    void toDefault();
  };

  class NALUnit
  {
    public:
      NALUnit(NALHeader header);
      virtual ~NALUnit();
      virtual NALUnitType getType() const;

      std::shared_ptr<NALUnit> copy() const;

      bool            m_processFailed;

      NALHeader     m_nalHeader;
  };



  class VPS: public NALUnit
  {
    public:
      VPS();
      uint8_t                   vps_video_parameter_set_id;
      uint8_t                   vps_base_layer_internal_flag;
      uint8_t                   vps_base_layer_available_flag;

      uint8_t                   vps_max_layers_minus1;
      uint8_t                   vps_max_sub_layers_minus1;
      uint8_t                   vps_temporal_id_nesting_flag;
      ProfileTierLevel          profile_tier_level;
      uint8_t                   vps_sub_layer_ordering_info_present_flag;
      std::vector<uint32_t>     vps_max_dec_pic_buffering_minus1;
      std::vector<uint32_t>     vps_max_num_reorder_pics;
      std::vector<uint32_t>     vps_max_latency_increase_plus1;
      uint8_t                   vps_max_layer_id;
      uint32_t                  vps_num_layer_sets_minus1;
      std::vector<std::vector<uint8_t> >
                                layer_id_included_flag;
      uint8_t                   vps_timing_info_present_flag;
      uint32_t                  vps_num_units_in_tick;
      uint32_t                  vps_time_scale;
      uint8_t                   vps_poc_proportional_to_timing_flag;
      uint32_t                  vps_num_ticks_poc_diff_one_minus1;
      uint32_t                  vps_num_hrd_parameters;
      std::vector<uint32_t>     hrd_layer_set_idx;
      std::vector<uint8_t>      cprms_present_flag;
      std::vector<HrdParameters>
                                hrd_parameters;
      uint8_t                   vps_extension_flag;

      VpsExtension vps_extension;

      uint8_t vps_extension2_flag;
      uint8_t vps_3d_extension_flag;
      Vps3dExtension vps_3d_extension;
      
      uint8_t vps_extension3_flag;

      void toDefault();
      bool operator == (const VPS &) const;
  };


  class SPS: public NALUnit
  {
    public:
      SPS();
      uint8_t                  sps_video_parameter_set_id;
      uint8_t                  sps_max_sub_layers_minus1;
      uint8_t                  sps_temporal_id_nesting_flag;
      ProfileTierLevel         profile_tier_level;
      uint32_t                 sps_seq_parameter_set_id;
      uint32_t                 chroma_format_idc;
      uint8_t                  separate_colour_plane_flag;
      uint32_t                 pic_width_in_luma_samples;
      uint32_t                 pic_height_in_luma_samples;
      uint8_t                  conformance_window_flag;
      uint32_t                 conf_win_left_offset;
      uint32_t                 conf_win_right_offset;
      uint32_t                 conf_win_top_offset;
      uint32_t                 conf_win_bottom_offset;
      uint32_t                 bit_depth_luma_minus8;
      uint32_t                 bit_depth_chroma_minus8;
      uint32_t                 log2_max_pic_order_cnt_lsb_minus4;
      uint8_t                  sps_sub_layer_ordering_info_present_flag;
      std::vector<uint32_t>    sps_max_dec_pic_buffering_minus1;
      std::vector<uint32_t>    sps_max_num_reorder_pics;
      std::vector<uint32_t>    sps_max_latency_increase_plus1;
      uint32_t                 log2_min_luma_coding_block_size_minus3;
      uint32_t                 log2_diff_max_min_luma_coding_block_size;
      uint32_t                 log2_min_transform_block_size_minus2;
      uint32_t                 log2_diff_max_min_transform_block_size;
      uint32_t                 max_transform_hierarchy_depth_inter;
      uint32_t                 max_transform_hierarchy_depth_intra;
      uint8_t                  scaling_list_enabled_flag;
      uint8_t                  sps_scaling_list_data_present_flag;
      ScalingListData          scaling_list_data;
      uint8_t                  amp_enabled_flag;
      uint8_t                  sample_adaptive_offset_enabled_flag;
      uint8_t                  pcm_enabled_flag;
      uint8_t                  pcm_sample_bit_depth_luma_minus1;
      uint8_t                  pcm_sample_bit_depth_chroma_minus1;
      uint32_t                 log2_min_pcm_luma_coding_block_size_minus3;
      uint32_t                 log2_diff_max_min_pcm_luma_coding_block_size;
      uint8_t                  pcm_loop_filter_disabled_flag;
      uint32_t                 num_short_term_ref_pic_sets;
      std::vector<ShortTermRefPicSet>
                               short_term_ref_pic_set;
      uint8_t                  long_term_ref_pics_present_flag;
      uint32_t                 num_long_term_ref_pics_sps;
      std::vector<uint32_t>    lt_ref_pic_poc_lsb_sps;
      std::vector<uint8_t>     used_by_curr_pic_lt_sps_flag;
      uint8_t                  sps_temporal_mvp_enabled_flag;
      uint8_t                  strong_intra_smoothing_enabled_flag;
      uint8_t                  vui_parameters_present_flag;
      VuiParameters            vui_parameters;
      uint8_t                  sps_extension_flag;

      uint8_t                  sps_range_extension_flag;
      uint8_t                  sps_multilayer_extension_flag;
      uint8_t                  sps_3d_extension_flag;
      uint8_t                  sps_extension_5bits;

      SpsRangeExtension sps_range_extension;
      uint8_t                  inter_view_mv_vert_constraint_flag;
      Sps3dExtension sps_3d_extension;

      void toDefault();

      bool operator == (const SPS &) const;
  };


  class PPS: public NALUnit
  {
    public:
    PPS();

    uint32_t     pps_pic_parameter_set_id;
    uint32_t     pps_seq_parameter_set_id;
    uint8_t      dependent_slice_segments_enabled_flag;
    uint8_t      output_flag_present_flag;
    uint8_t      num_extra_slice_header_bits;
    uint8_t      sign_data_hiding_flag;
    uint8_t      cabac_init_present_flag;
    uint32_t     num_ref_idx_l0_default_active_minus1;
    uint32_t     num_ref_idx_l1_default_active_minus1;
    int32_t      init_qp_minus26;
    uint8_t      constrained_intra_pred_flag;
    uint8_t      transform_skip_enabled_flag;
    uint8_t      cu_qp_delta_enabled_flag;
    uint32_t     diff_cu_qp_delta_depth;
    int32_t      pps_cb_qp_offset;
    int32_t      pps_cr_qp_offset;
    uint8_t      pps_slice_chroma_qp_offsets_present_flag;
    uint8_t      weighted_pred_flag;
    uint8_t      weighted_bipred_flag;
    uint8_t      transquant_bypass_enabled_flag;
    uint8_t      tiles_enabled_flag;
    uint8_t      entropy_coding_sync_enabled_flag;
    uint32_t     num_tile_columns_minus1;
    uint32_t     num_tile_rows_minus1;
    uint8_t      uniform_spacing_flag;
    std::vector<uint32_t>
                 column_width_minus1;
    std::vector<uint32_t>
                 row_height_minus1;
    uint8_t      loop_filter_across_tiles_enabled_flag;
    uint8_t      pps_loop_filter_across_slices_enabled_flag;
    uint8_t      deblocking_filter_control_present_flag;
    uint8_t      deblocking_filter_override_enabled_flag;
    uint8_t      pps_deblocking_filter_disabled_flag;
    uint32_t     pps_beta_offset_div2;
    uint32_t     pps_tc_offset_div2;
    uint8_t      pps_scaling_list_data_present_flag;
    ScalingListData
                 scaling_list_data;
    uint8_t      lists_modification_present_flag;
    int32_t      log2_parallel_merge_level_minus2;
    uint8_t      slice_segment_header_extension_present_flag;
    uint8_t      pps_extension_flag;

    uint8_t pps_range_extension_flag;
    uint8_t pps_multilayer_extension_flag;
    uint8_t pps_3d_extension_flag;
    uint8_t pps_extension_5bits;

    PpsRangeExtension pps_range_extension;

    PpsMultilayerExtension pps_multilayer_extension;

    Pps3dExtension pps_3d_extension;

    void toDefault();

    bool operator == (const PPS &) const;
  };



  class Slice: public NALUnit
  {
    public:
      enum SliceType
      {
        B_SLICE = 0,
        P_SLICE = 1,
        I_SLICE = 2,
        NONE_SLICE = 3
      };

      Slice(NALHeader header);
      uint8_t                  first_slice_segment_in_pic_flag;
      uint8_t                  no_output_of_prior_pics_flag;
      uint32_t                 slice_pic_parameter_set_id;
      uint8_t                  dependent_slice_segment_flag;
      uint32_t                 slice_segment_address;
      std::vector<uint32_t>    slice_reserved_undetermined_flag;
      uint32_t                 slice_type;
      uint8_t                  pic_output_flag;
      uint8_t                  colour_plane_id;
      uint32_t                 slice_pic_order_cnt_lsb;
      uint8_t                  short_term_ref_pic_set_sps_flag;
      ShortTermRefPicSet       short_term_ref_pic_set;
      uint8_t                  short_term_ref_pic_set_idx;

      uint32_t                 num_long_term_sps;
      uint32_t                 num_long_term_pics;
      std::vector<uint32_t>    lt_idx_sps;
      std::vector<uint32_t>    poc_lsb_lt;
      std::vector<uint8_t>     used_by_curr_pic_lt_flag;
      std::vector<uint8_t>     delta_poc_msb_present_flag;
      std::vector<uint32_t>    delta_poc_msb_cycle_lt;

      uint8_t                  slice_temporal_mvp_enabled_flag;
      uint8_t                  slice_sao_luma_flag;
      uint8_t                  slice_sao_chroma_flag;
      uint8_t                  num_ref_idx_active_override_flag;
      uint32_t                 num_ref_idx_l0_active_minus1;
      uint32_t                 num_ref_idx_l1_active_minus1;
      RefPicListModification   ref_pic_lists_modification;
      uint8_t                  mvd_l1_zero_flag;
      uint8_t                  cabac_init_flag;
      uint8_t                  collocated_from_l0_flag;
      uint32_t                 collocated_ref_idx;
      PredWeightTable          pred_weight_table;
      uint32_t                 five_minus_max_num_merge_cand;
      int32_t                  slice_qp_delta;
      int32_t                  slice_cb_qp_offset;
      int32_t                  slice_cr_qp_offset;
      uint8_t                  deblocking_filter_override_flag;
      uint8_t                  slice_deblocking_filter_disabled_flag;
      int32_t                  slice_beta_offset_div2;
      int32_t                  slice_tc_offset_div2;
      int32_t                  slice_loop_filter_across_slices_enabled_flag;
      uint32_t                 num_entry_point_offsets;
      uint32_t                 offset_len_minus1;
      std::vector<uint32_t>    entry_point_offset_minus1;
      uint32_t                 slice_segment_header_extension_length;
      std::vector<uint8_t>     slice_segment_header_extension_data_byte;


      void toDefault();
  };

  class AUD: public NALUnit
  {
  public:
    AUD();

    uint8_t            pic_type;
    void toDefault();
  };


  class SEI: public NALUnit
  {
  public:
    SEI(NALHeader header);
    std::vector<SeiMessage>     sei_message;

    void toDefault();
  };

}

#endif
