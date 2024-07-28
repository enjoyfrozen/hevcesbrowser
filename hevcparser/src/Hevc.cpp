#include "Hevc.h"

#include <string.h>
#include <iostream>
#include <thread>


using namespace HEVC;

NALUnit::NALUnit(NALHeader header):
  m_nalHeader(header)
  ,m_processFailed(false)
{
    //std::cout << __FUNCTION__ << " nal type:" << m_nalHeader.type << " pointer "<<this <<" thread "<< std::this_thread::get_id() << std::endl;
}


NALUnit::~NALUnit()
{
    //std::cout << __FUNCTION__ << " nal type:"<< m_nalHeader.type << " pointer " << this << " thread " << std::this_thread::get_id() <<std::endl;
}


NALUnitType NALUnit::getType() const
{
  return m_nalHeader.type;
}


std::shared_ptr<NALUnit> NALUnit::copy() const
{
  std::shared_ptr<NALUnit> res;

  switch(m_nalHeader.type)
  {
    case NAL_VPS:
    {
      res = std::shared_ptr<NALUnit>(new VPS(*((VPS *) this)));
      break;
    }

    case NAL_SPS:
    {
      res = std::shared_ptr<NALUnit>(new SPS(*((SPS *) this)));
      break;
    }

    case NAL_PPS:
    {
      res = std::shared_ptr<NALUnit>(new PPS(*((PPS *) this)));
      break;
    }

    case NAL_TRAIL_R:
    case NAL_TRAIL_N:
    case NAL_TSA_N:
    case NAL_TSA_R:
    case NAL_STSA_N:
    case NAL_STSA_R:
    case NAL_BLA_W_LP:
    case NAL_BLA_W_RADL:
    case NAL_BLA_N_LP:
    case NAL_IDR_W_RADL:
    case NAL_IDR_N_LP:
    case NAL_CRA_NUT:
    case NAL_RADL_N:
    case NAL_RADL_R:
    case NAL_RASL_N:
    case NAL_RASL_R:
    {
      res = std::shared_ptr<NALUnit>(new Slice(*((Slice *) this)));
      break;
    }

    case NAL_SEI_PREFIX:
    {
        res = std::shared_ptr<NALUnit>(new SEI(*((SEI*)this)));
        break;
    }


    default:
      res = std::shared_ptr<NALUnit>(new NALUnit(*this));
  };

  res -> m_processFailed = m_processFailed;
  return res;
}

HEVC::VPS::VPS(): NALUnit({HEVC::NAL_VPS, 0, 0})
{
  toDefault();
}


HEVC::SPS::SPS(): NALUnit({HEVC::NAL_SPS, 0, 0})
{
  toDefault();
}


HEVC::PPS::PPS(): NALUnit({HEVC::NAL_PPS, 0, 0})
{
  toDefault();
};


HEVC::AUD::AUD(): NALUnit({HEVC::NAL_AUD, 0, 0})
{
  toDefault();
};


HEVC::SEI::SEI(NALHeader header): NALUnit(header)
{
  toDefault();
};


HEVC::Slice::Slice(NALHeader header): NALUnit(header)
{
  toDefault();
};


bool HEVC::VPS::operator == (const HEVC::VPS &vps) const
{
  if(this == &vps)
    return true;

  if(vps_video_parameter_set_id != vps.vps_video_parameter_set_id)
    return false;

  if(vps_max_layers_minus1 != vps.vps_max_layers_minus1)
    return false;

  if(vps_max_sub_layers_minus1 != vps.vps_max_sub_layers_minus1)
    return false;

  if(vps_temporal_id_nesting_flag != vps.vps_temporal_id_nesting_flag)
    return false;

  if(!(profile_tier_level == vps.profile_tier_level))
    return false;

  if(vps_sub_layer_ordering_info_present_flag != vps.vps_sub_layer_ordering_info_present_flag)
    return false;

  if(vps_max_dec_pic_buffering_minus1 != vps.vps_max_dec_pic_buffering_minus1)
    return false;

  if(vps_max_num_reorder_pics != vps.vps_max_num_reorder_pics)
    return false;

  if(vps_max_latency_increase_plus1 != vps.vps_max_latency_increase_plus1)
    return false;

  if(vps_max_layer_id != vps.vps_max_layer_id)
    return false;

  if(vps_num_layer_sets_minus1 != vps.vps_num_layer_sets_minus1)
    return false;

  if(layer_id_included_flag != vps.layer_id_included_flag)
    return false;

  if(vps_timing_info_present_flag != vps.vps_timing_info_present_flag)
    return false;

  if(vps_num_units_in_tick != vps.vps_num_units_in_tick)
    return false;

  if(vps_time_scale != vps.vps_time_scale)
    return false;

  if(vps_poc_proportional_to_timing_flag != vps.vps_poc_proportional_to_timing_flag)
    return false;

  if(vps_num_ticks_poc_diff_one_minus1 != vps.vps_num_ticks_poc_diff_one_minus1)
    return false;

  if(vps_num_hrd_parameters != vps.vps_num_hrd_parameters)
    return false;

  if(hrd_layer_set_idx!= vps.hrd_layer_set_idx)
    return false;

  if(cprms_present_flag != vps.cprms_present_flag)
    return false;

  if(!(hrd_parameters == vps.hrd_parameters))
    return false;

  if(vps_extension_flag != vps.vps_extension_flag)
    return false;

  return true;
}


bool HEVC::SPS::operator == (const HEVC::SPS &sps) const
{
  if(this == &sps)
    return true;

  if(sps_video_parameter_set_id != sps.sps_video_parameter_set_id)
    return false;

  if(sps_max_sub_layers_minus1 != sps.sps_max_sub_layers_minus1)
    return false;

  if(sps_temporal_id_nesting_flag != sps.sps_temporal_id_nesting_flag)
    return false;

  if(!(profile_tier_level == sps.profile_tier_level))
    return false;

  if(sps_seq_parameter_set_id != sps.sps_seq_parameter_set_id)
    return false;

  if(chroma_format_idc != sps.chroma_format_idc)
    return false;

  if(separate_colour_plane_flag != sps.separate_colour_plane_flag)
    return false;

  if(pic_width_in_luma_samples != sps.pic_width_in_luma_samples)
    return false;

  if(pic_height_in_luma_samples != sps.pic_height_in_luma_samples)
    return false;

  if(conformance_window_flag != sps.conformance_window_flag)
    return false;

  if(conf_win_left_offset != sps.conf_win_left_offset)
    return false;

  if(conf_win_right_offset != sps.conf_win_right_offset)
    return false;

  if(conf_win_top_offset != sps.conf_win_top_offset)
    return false;

  if(conf_win_bottom_offset != sps.conf_win_bottom_offset)
    return false;

  if(bit_depth_luma_minus8 != sps.bit_depth_luma_minus8)
    return false;

  if(bit_depth_chroma_minus8 != sps.bit_depth_chroma_minus8)
    return false;

  if(log2_max_pic_order_cnt_lsb_minus4 != sps.log2_max_pic_order_cnt_lsb_minus4)
    return false;

  if(sps_sub_layer_ordering_info_present_flag != sps.sps_sub_layer_ordering_info_present_flag)
    return false;

  if(sps_max_dec_pic_buffering_minus1 != sps.sps_max_dec_pic_buffering_minus1)
    return false;

  if(sps_max_num_reorder_pics != sps.sps_max_num_reorder_pics)
    return false;

  if(sps_max_latency_increase_plus1 != sps.sps_max_latency_increase_plus1)
    return false;

  if(log2_min_luma_coding_block_size_minus3 != sps.log2_min_luma_coding_block_size_minus3)
    return false;

  if(log2_diff_max_min_luma_coding_block_size != sps.log2_diff_max_min_luma_coding_block_size)
    return false;

  if(log2_min_transform_block_size_minus2 != sps.log2_min_transform_block_size_minus2)
    return false;

  if(log2_diff_max_min_transform_block_size != sps.log2_diff_max_min_transform_block_size)
    return false;

  if(max_transform_hierarchy_depth_inter != sps.max_transform_hierarchy_depth_inter)
    return false;

  if(max_transform_hierarchy_depth_intra != sps.max_transform_hierarchy_depth_intra)
    return false;

  if(scaling_list_enabled_flag != sps.scaling_list_enabled_flag)
    return false;

  if(!(scaling_list_data == sps.scaling_list_data))
    return false;

  if(sps_scaling_list_data_present_flag != sps.sps_scaling_list_data_present_flag)
    return false;

  if(amp_enabled_flag != sps.amp_enabled_flag)
    return false;

  if(sample_adaptive_offset_enabled_flag != sps.sample_adaptive_offset_enabled_flag)
    return false;

  if(pcm_enabled_flag != sps.pcm_enabled_flag)
    return false;

  if(num_short_term_ref_pic_sets != sps.num_short_term_ref_pic_sets)
    return false;

  if(short_term_ref_pic_set != sps.short_term_ref_pic_set)
    return false;

  if(long_term_ref_pics_present_flag != sps.long_term_ref_pics_present_flag)
    return false;

  if(num_long_term_ref_pics_sps != sps.num_long_term_ref_pics_sps)
    return false;

  if(lt_ref_pic_poc_lsb_sps != sps.lt_ref_pic_poc_lsb_sps)
    return false;

  if(used_by_curr_pic_lt_sps_flag != sps.used_by_curr_pic_lt_sps_flag)
    return false;

  if(sps_temporal_mvp_enabled_flag != sps.sps_temporal_mvp_enabled_flag)
    return false;

  if(strong_intra_smoothing_enabled_flag != sps.strong_intra_smoothing_enabled_flag)
    return false;

  if(vui_parameters_present_flag != sps.vui_parameters_present_flag)
    return false;

  if(!(vui_parameters == sps.vui_parameters))
    return false;

  if(sps_extension_flag != sps.sps_extension_flag)
    return false;

  return true;
}


bool HEVC::PPS::operator == (const HEVC::PPS &pps) const
{
  if(this == &pps)
    return true;

  if(pps_pic_parameter_set_id != pps.pps_pic_parameter_set_id)
    return false;

  if(pps_seq_parameter_set_id != pps.pps_seq_parameter_set_id)
    return false;

  if(dependent_slice_segments_enabled_flag != pps.dependent_slice_segments_enabled_flag)
    return false;

  if(output_flag_present_flag != pps.output_flag_present_flag)
    return false;

  if(num_extra_slice_header_bits != pps.num_extra_slice_header_bits)
    return false;

  if(sign_data_hiding_flag != pps.sign_data_hiding_flag)
    return false;

  if(cabac_init_present_flag != pps.cabac_init_present_flag)
    return false;

  if(num_ref_idx_l0_default_active_minus1 != pps.num_ref_idx_l0_default_active_minus1)
    return false;

  if(num_ref_idx_l1_default_active_minus1 != pps.num_ref_idx_l1_default_active_minus1)
    return false;

  if(init_qp_minus26 != pps.init_qp_minus26)
    return false;

  if(constrained_intra_pred_flag != pps.constrained_intra_pred_flag)
    return false;

  if(transform_skip_enabled_flag != pps.transform_skip_enabled_flag)
    return false;

  if(cu_qp_delta_enabled_flag != pps.cu_qp_delta_enabled_flag)
    return false;

  if(diff_cu_qp_delta_depth != pps.diff_cu_qp_delta_depth)
    return false;

  if(pps_cb_qp_offset != pps.pps_cb_qp_offset)
    return false;

  if(pps_cr_qp_offset != pps.pps_cr_qp_offset)
    return false;

  if(pps_slice_chroma_qp_offsets_present_flag != pps.pps_slice_chroma_qp_offsets_present_flag)
    return false;

  if(weighted_pred_flag != pps.weighted_pred_flag)
    return false;

  if(weighted_bipred_flag != pps.weighted_bipred_flag)
    return false;

  if(transquant_bypass_enabled_flag != pps.transquant_bypass_enabled_flag)
    return false;

  if(tiles_enabled_flag != pps.tiles_enabled_flag)
    return false;

  if(entropy_coding_sync_enabled_flag != pps.entropy_coding_sync_enabled_flag)
    return false;

  if(num_tile_columns_minus1 != pps.num_tile_columns_minus1)
    return false;

  if(num_tile_rows_minus1 != pps.num_tile_rows_minus1)
    return false;

  if(uniform_spacing_flag != pps.uniform_spacing_flag)
    return false;

  if(column_width_minus1 != pps.column_width_minus1)
    return false;

  if(row_height_minus1 != pps.row_height_minus1)
    return false;

  if(loop_filter_across_tiles_enabled_flag != pps.loop_filter_across_tiles_enabled_flag)
    return false;

  if(pps_loop_filter_across_slices_enabled_flag != pps.pps_loop_filter_across_slices_enabled_flag)
    return false;

  if(deblocking_filter_control_present_flag != pps.deblocking_filter_control_present_flag)
    return false;

  if(deblocking_filter_override_enabled_flag != pps.deblocking_filter_override_enabled_flag)
    return false;

  if(pps_deblocking_filter_disabled_flag != pps.pps_deblocking_filter_disabled_flag)
    return false;

  if(pps_beta_offset_div2 != pps.pps_beta_offset_div2)
    return false;

  if(pps_tc_offset_div2 != pps.pps_tc_offset_div2)
    return false;

  if(pps_scaling_list_data_present_flag != pps_scaling_list_data_present_flag)
    return false;

  if(lists_modification_present_flag != lists_modification_present_flag)
    return false;

  if(log2_parallel_merge_level_minus2 != log2_parallel_merge_level_minus2)
    return false;

  if(slice_segment_header_extension_present_flag != slice_segment_header_extension_present_flag)
    return false;

  if(pps_extension_flag != pps_extension_flag)
    return false;

  return true;
}



bool HEVC::ProfileTierLevel::operator == (const HEVC::ProfileTierLevel &ptl) const
{
  if(this == &ptl)
    return true;

  if(general_profile_space != ptl.general_profile_space)
    return false;

  if(general_tier_flag != ptl.general_tier_flag)
    return false;

  if(general_profile_idc != ptl.general_profile_idc)
    return false;

  if(general_progressive_source_flag != ptl.general_progressive_source_flag)
    return false;

  if(general_interlaced_source_flag != ptl.general_interlaced_source_flag)
    return false;

  if(general_non_packed_constraint_flag != ptl.general_non_packed_constraint_flag)
    return false;

  if(general_frame_only_constraint_flag != ptl.general_frame_only_constraint_flag)
    return false;

  if(general_level_idc != ptl.general_level_idc)
    return false;

  if(sub_layer_profile_present_flag != ptl.sub_layer_profile_present_flag)
    return false;

  if(sub_layer_level_present_flag != ptl.sub_layer_level_present_flag)
    return false;

  if(memcmp(general_profile_compatibility_flag, ptl.general_profile_compatibility_flag, 32) != 0)
    return false;

  if(sub_layer_profile_space != ptl.sub_layer_profile_space)
    return false;

  if(sub_layer_tier_flag != ptl.sub_layer_tier_flag)
    return false;

  if(sub_layer_profile_idc != ptl.sub_layer_profile_idc)
    return false;

  if(sub_layer_profile_compatibility_flag != ptl.sub_layer_profile_compatibility_flag)
    return false;

  if(sub_layer_progressive_source_flag != ptl.sub_layer_progressive_source_flag)
    return false;

  if(sub_layer_interlaced_source_flag != ptl.sub_layer_interlaced_source_flag)
    return false;

  if(sub_layer_non_packed_constraint_flag != ptl.sub_layer_non_packed_constraint_flag)
    return false;

  if(sub_layer_frame_only_constraint_flag != ptl.sub_layer_frame_only_constraint_flag)
    return false;

  if(sub_layer_level_idc != ptl.sub_layer_level_idc)
    return false;

  return true;
}



bool HEVC::HrdParameters::operator == (const HrdParameters &hrd) const
{
  if(this == &hrd)
    return true;

  if(nal_hrd_parameters_present_flag != hrd.nal_hrd_parameters_present_flag)
    return false;

  if(vcl_hrd_parameters_present_flag != hrd.vcl_hrd_parameters_present_flag)
    return false;

  if(sub_pic_hrd_params_present_flag != hrd.sub_pic_hrd_params_present_flag)
    return false;

  if(tick_divisor_minus2 != hrd.tick_divisor_minus2)
    return false;

  if(du_cpb_removal_delay_increment_length_minus1 != hrd.du_cpb_removal_delay_increment_length_minus1)
    return false;

  if(sub_pic_cpb_params_in_pic_timing_sei_flag != hrd.sub_pic_cpb_params_in_pic_timing_sei_flag)
    return false;

  if(dpb_output_delay_du_length_minus1 != hrd.dpb_output_delay_du_length_minus1)
    return false;

  if(bit_rate_scale != hrd.bit_rate_scale)
    return false;

  if(cpb_size_scale != hrd.cpb_size_scale)
    return false;

  if(cpb_size_du_scale != hrd.cpb_size_du_scale)
    return false;

  if(initial_cpb_removal_delay_length_minus1 != hrd.initial_cpb_removal_delay_length_minus1)
    return false;

  if(au_cpb_removal_delay_length_minus1 != hrd.au_cpb_removal_delay_length_minus1)
    return false;

  if(dpb_output_delay_length_minus1 != hrd.dpb_output_delay_length_minus1)
    return false;

  if(fixed_pic_rate_general_flag != hrd.fixed_pic_rate_general_flag)
    return false;

  if(fixed_pic_rate_within_cvs_flag != hrd.fixed_pic_rate_within_cvs_flag)
    return false;

  if(elemental_duration_in_tc_minus1 != hrd.elemental_duration_in_tc_minus1)
    return false;

  if(low_delay_hrd_flag != hrd.low_delay_hrd_flag)
    return false;

  if(cpb_cnt_minus1 != hrd.cpb_cnt_minus1)
    return false;

  if(!(nal_sub_layer_hrd_parameters == hrd.nal_sub_layer_hrd_parameters))
    return false;

  if(!(vcl_sub_layer_hrd_parameters == hrd.vcl_sub_layer_hrd_parameters))
    return false;

  return true;
}



bool HEVC::SubLayerHrdParameters::operator == (const SubLayerHrdParameters &slhrd) const
{
  if(this == &slhrd)
    return true;

  if(bit_rate_value_minus1 != slhrd.bit_rate_value_minus1)
    return false;

  if(cpb_size_value_minus1 != slhrd.cpb_size_value_minus1)
    return false;

  if(cpb_size_du_value_minus1 != slhrd.cpb_size_du_value_minus1)
    return false;

  if(bit_rate_du_value_minus1 != slhrd.bit_rate_du_value_minus1)
    return false;

  if(cbr_flag != slhrd.cbr_flag)
    return false;

  return true;
}


bool HEVC::ShortTermRefPicSet::operator == (const ShortTermRefPicSet &strps) const
{
  if(this == &strps)
    return true;

  if(inter_ref_pic_set_prediction_flag != strps.inter_ref_pic_set_prediction_flag)
    return false;

  if(delta_idx_minus1 != strps.delta_idx_minus1)
    return false;

  if(delta_rps_sign != strps.delta_rps_sign)
    return false;

  if(abs_delta_rps_minus1 != strps.abs_delta_rps_minus1)
    return false;

  if(used_by_curr_pic_flag != strps.used_by_curr_pic_flag)
    return false;

  if(use_delta_flag != strps.use_delta_flag)
    return false;

  if(num_negative_pics != strps.num_negative_pics)
    return false;

  if(num_positive_pics != strps.num_positive_pics)
    return false;

  if(delta_poc_s0_minus1 != strps.delta_poc_s0_minus1)
    return false;

  if(used_by_curr_pic_s0_flag != strps.used_by_curr_pic_s0_flag)
    return false;

  if(delta_poc_s1_minus1 != strps.delta_poc_s1_minus1)
    return false;

  if(used_by_curr_pic_s1_flag != strps.used_by_curr_pic_s1_flag)
    return false;

  return true;
}



bool HEVC::VuiParameters::operator == (const VuiParameters &vui) const
{
  if(this == &vui)
    return true;

  if(aspect_ratio_info_present_flag != vui.aspect_ratio_info_present_flag)
    return false;

  if(aspect_ratio_idc != vui.aspect_ratio_idc)
    return false;

  if(sar_width != vui.sar_width)
    return false;

  if(sar_height != vui.sar_height)
    return false;

  if(overscan_info_present_flag != vui.overscan_info_present_flag)
    return false;

  if(overscan_appropriate_flag != vui.overscan_appropriate_flag)
    return false;

  if(video_signal_type_present_flag != vui.video_signal_type_present_flag)
    return false;

  if(video_format != vui.video_format)
    return false;

  if(video_full_range_flag != vui.video_full_range_flag)
    return false;

  if(colour_description_present_flag != vui.colour_description_present_flag)
    return false;

  if(colour_primaries != vui.colour_primaries)
    return false;

  if(transfer_characteristics != vui.transfer_characteristics)
    return false;

  if(matrix_coeffs != vui.matrix_coeffs)
    return false;

  if(chroma_loc_info_present_flag != vui.chroma_loc_info_present_flag)
    return false;

  if(chroma_sample_loc_type_top_field != vui.chroma_sample_loc_type_top_field)
    return false;

  if(chroma_sample_loc_type_bottom_field != vui.chroma_sample_loc_type_bottom_field)
    return false;

  if(neutral_chroma_indication_flag != vui.neutral_chroma_indication_flag)
    return false;

  if(field_seq_flag != vui.field_seq_flag)
    return false;

  if(frame_field_info_present_flag != vui.frame_field_info_present_flag)
    return false;

  if(default_display_window_flag != vui.default_display_window_flag)
    return false;

  if(def_disp_win_left_offset != vui.def_disp_win_left_offset)
    return false;

  if(def_disp_win_right_offset != vui.def_disp_win_right_offset)
    return false;

  if(def_disp_win_top_offset != vui.def_disp_win_top_offset)
    return false;

  if(def_disp_win_bottom_offset != vui.def_disp_win_bottom_offset)
    return false;

  if(vui_timing_info_present_flag != vui.vui_timing_info_present_flag)
    return false;

  if(vui_num_units_in_tick != vui.vui_num_units_in_tick)
    return false;

  if(vui_time_scale != vui.vui_time_scale)
    return false;

  if(vui_poc_proportional_to_timing_flag != vui.vui_poc_proportional_to_timing_flag)
    return false;

  if(vui_num_ticks_poc_diff_one_minus1 != vui.vui_num_ticks_poc_diff_one_minus1)
    return false;

  if(vui_hrd_parameters_present_flag != vui.vui_hrd_parameters_present_flag)
    return false;

  if(!(hrd_parameters == vui.hrd_parameters))
    return false;

  if(bitstream_restriction_flag != vui.bitstream_restriction_flag)
    return false;

  if(tiles_fixed_structure_flag != vui.tiles_fixed_structure_flag)
    return false;

  if(motion_vectors_over_pic_boundaries_flag != vui.motion_vectors_over_pic_boundaries_flag)
    return false;

  if(restricted_ref_pic_lists_flag != vui.restricted_ref_pic_lists_flag)
    return false;

  if(min_spatial_segmentation_idc != vui.min_spatial_segmentation_idc)
    return false;

  if(max_bytes_per_pic_denom != vui.max_bytes_per_pic_denom)
    return false;

  if(max_bits_per_min_cu_denom != vui.max_bits_per_min_cu_denom)
    return false;

  if(log2_max_mv_length_horizontal != vui.log2_max_mv_length_horizontal)
    return false;

  if(log2_max_mv_length_vertical != vui.log2_max_mv_length_vertical)
    return false;

  return true;
}


bool ScalingListData::operator == (const ScalingListData &obj) const
{
  if(this == &obj)
    return true;

  if(scaling_list_pred_mode_flag != obj.scaling_list_pred_mode_flag)
    return false;

  if(scaling_list_pred_matrix_id_delta != obj.scaling_list_pred_matrix_id_delta)
    return false;

  if(scaling_list_dc_coef_minus8 != obj.scaling_list_dc_coef_minus8)
    return false;

    if(scaling_list_delta_coef != obj.scaling_list_delta_coef)
    return false;

  return true;
}


bool RefPicListModification::operator == (const RefPicListModification &obj) const
{
  if(this == &obj)
    return true;

  if(ref_pic_list_modification_flag_l0 != ref_pic_list_modification_flag_l0)
    return false;

  if(list_entry_l0 != obj.list_entry_l0)
    return false;

  if(ref_pic_list_modification_flag_l1 != obj.ref_pic_list_modification_flag_l1)
    return false;

  if(list_entry_l1 != obj.list_entry_l1)
    return false;

  return true;
}


void ProfileTierLevel::toDefault()
{
    general_profile_space = 0;
    general_tier_flag = 0;
    general_profile_idc = 0;
    memset(&general_profile_compatibility_flag[0], 0, sizeof(general_profile_compatibility_flag));
    general_progressive_source_flag = 0;
    general_interlaced_source_flag = 0;
    general_non_packed_constraint_flag = 0;
    general_frame_only_constraint_flag = 0;

    general_max_12bit_constraint_flag = 0;
    general_max_10bit_constraint_flag = 0;
    general_max_8bit_constraint_flag = 0;
    general_max_422chroma_constraint_flag = 0;
    general_max_420chroma_constraint_flag = 0;
    general_max_monochrome_constraint_flag = 0;
    general_intra_constraint_flag = 0;
    general_one_picture_only_constraint_flag = 0;
    general_lower_bit_rate_constraint_flag = 0;
    general_inbld_flag = 0;


    general_level_idc = 0;
    sub_layer_profile_present_flag.clear();
    sub_layer_level_present_flag.clear();
    sub_layer_profile_space.clear();
    sub_layer_tier_flag.clear();
    sub_layer_profile_idc.clear();
    sub_layer_profile_compatibility_flag.clear();
    sub_layer_progressive_source_flag.clear();
    sub_layer_interlaced_source_flag.clear();
    sub_layer_non_packed_constraint_flag.clear();
    sub_layer_frame_only_constraint_flag.clear();
    sub_layer_level_idc.clear();

    sub_layer_max_12bit_constraint_flag.clear();
    sub_layer_max_10bit_constraint_flag.clear();
    sub_layer_max_8bit_constraint_flag.clear();
    sub_layer_max_422chroma_constraint_flag.clear();
    sub_layer_max_420chroma_constraint_flag.clear();
    sub_layer_max_monochrome_constraint_flag.clear();
    sub_layer_intra_constraint_flag.clear();
    sub_layer_one_picture_only_constraint_flag.clear();
    sub_layer_lower_bit_rate_constraint_flag.clear();
    sub_layer_inbld_flag.clear();
}


void SubLayerHrdParameters::toDefault()
{
    bit_rate_value_minus1.clear();
    cpb_size_value_minus1.clear();
    cpb_size_du_value_minus1.clear();
    bit_rate_du_value_minus1.clear();
    cbr_flag.clear();
}

void VideoSignalInfo::toDefault()
{
    video_vps_format = 0;
    video_full_range_vps_flag = 0;
    colour_primaries_vps = 0;
    transfer_characteristics_vps = 0;
    matrix_coeffs_vps = 0;
}

void HrdParameters::toDefault()
{
  nal_hrd_parameters_present_flag = 0;
  vcl_hrd_parameters_present_flag = 0;
  sub_pic_hrd_params_present_flag = 0;
  tick_divisor_minus2 = 0;
  du_cpb_removal_delay_increment_length_minus1 = 0;
  sub_pic_cpb_params_in_pic_timing_sei_flag = 0;
  dpb_output_delay_du_length_minus1 = 0;
  bit_rate_scale = 0;
  cpb_size_scale = 0;
  cpb_size_du_scale = 0;
  initial_cpb_removal_delay_length_minus1 = 23;
  au_cpb_removal_delay_length_minus1 = 23;
  dpb_output_delay_length_minus1 = 23;
  fixed_pic_rate_general_flag.clear();
  fixed_pic_rate_within_cvs_flag.clear();
  elemental_duration_in_tc_minus1.clear();
  low_delay_hrd_flag.clear();
  cpb_cnt_minus1.clear();
  nal_sub_layer_hrd_parameters.clear();
  vcl_sub_layer_hrd_parameters.clear();
}

void Sps3dExtension::toDefault()
{
    for (int i = 0; i <= 1;i++) {
         iv_di_mc_enabled_flag[i]=0;
         iv_mv_scal_enabled_flag[i] = 0;
         log2_ivmc_sub_pb_size_minus3[i] = 0;
         iv_res_pred_enabled_flag[i] = 0;
         depth_ref_enabled_flag[i] = 0;
         vsp_mc_enabled_flag[i] = 0;
         dbbp_enabled_flag[i] = 0;
         tex_mc_enabled_flag[i] = 0;
         log2_texmc_sub_pb_size_minus3[i] = 0;
         intra_contour_enabled_flag[i] = 0;
         intra_dc_only_wedge_enabled_flag[i] = 0;
         cqt_cu_part_pred_enabled_flag[i] = 0;
         inter_dc_only_enabled_flag[i] = 0;
         skip_intra_enabled_flag[i] = 0;
    }
}

void SpsRangeExtension::toDefault()
{
    transform_skip_rotation_enabled_flag = 0;
    transform_skip_context_enabled_flag = 0;
    implicit_rdpcm_enabled_flag = 0;
    explicit_rdpcm_enabled_flag = 0;

    extended_precision_processing_flag = 0;

    intra_smoothing_disabled_flag = 0;
    high_precision_offsets_enabled_flag = 0;
    persistent_rice_adaptation_enabled_flag = 0;
    cabac_bypass_alignment_enabled_flag = 0;
}

void PpsRangeExtension::toDefault()
{
    log2_max_transform_skip_block_size_minus2 = 0;
    cross_component_prediction_enabled_flag = 0;
    chroma_qp_offset_list_enabled_flag = 0;
    diff_cu_chroma_qp_offset_depth = 0;
    chroma_qp_offset_list_len_minus1 = 0;
    cb_qp_offset_list.clear();
    cr_qp_offset_list.clear();

    log2_sao_offset_scale_luma = 0;
    log2_sao_offset_scale_chroma = 0;
}

void ColourMappingTable::toDefault()
{
    num_cm_ref_layers_minus1 = 0;
    cm_ref_layer_id.clear();
    cm_octant_depth = 0;
    cm_y_part_num_log2 = 0;
    luma_bit_depth_cm_input_minus8 = 0;
    chroma_bit_depth_cm_input_minus8 = 0;
    luma_bit_depth_cm_output_minus8 = 0;
    chroma_bit_depth_cm_output_minus8 = 0;
    cm_res_quant_bits = 0;
    cm_delta_flc_bits_minus1 = 0;
    cm_adapt_threshold_u_delta = 0;
    cm_adapt_threshold_v_delta = 0;
}

void PpsMultilayerExtension::toDefault()
{
    poc_reset_info_present_flag = 0;
    pps_infer_scaling_list_flag = 0;
    pps_scaling_list_ref_layer_id = 0;
    num_ref_loc_offsets = 0;
    //for (i = 0; i < num_ref_loc_offsets; i++) 
    //{
    ref_loc_offset_layer_id.clear();
    scaled_ref_layer_offset_present_flag.clear();
    scaled_ref_layer_left_offset.clear();
    scaled_ref_layer_top_offset.clear();
    scaled_ref_layer_right_offset.clear();
    scaled_ref_layer_bottom_offset.clear();

    //}
    ref_region_offset_present_flag.clear();
    ref_region_left_offset.clear();
    ref_region_top_offset.clear();
    ref_region_right_offset.clear();
    ref_region_bottom_offset.clear();
    resample_phase_set_present_flag.clear();
    phase_hor_luma.clear();
    phase_ver_luma.clear();
    phase_hor_chroma_plus8.clear();
    phase_ver_chroma_plus8.clear();
    colour_mapping_enabled_flag = 0;
}

void Pps3dExtension::toDefault()
{
    dlts_present_flag = 0;
    pps_depth_layers_minus1 = 0;
    pps_bit_depth_for_depth_layers_minus8 = 0;
    dlt_flag.clear();
    dlt_pred_flag.clear();
    dlt_val_flags_present_flag.clear();
    dlt_value_flag.clear();
    delta_dlt.clear();
}


void DeltaDlt::toDefault()
{
    num_val_delta_dlt = 0;
    max_diff = 0;
    min_diff_minus1 = 0;
    delta_dlt_val0 = 0;
    delta_val_diff_minus_min.clear();
}


void Vps3dExtension::toDefault()
{
    cp_precision = 0;

    cp_in_slice_segment_header_flag.clear();
    cp_ref_voi.clear();
    vps_cp_scale.clear();
    vps_cp_off.clear();
    vps_cp_inv_scale_plus_scale.clear();
    vps_cp_inv_off_plus_off.clear();

}

void ShortTermRefPicSet::toDefault()
{
    inter_ref_pic_set_prediction_flag = 0;
    delta_idx_minus1 = 0;
    delta_rps_sign = 0;
    abs_delta_rps_minus1 = 0;
    used_by_curr_pic_flag.clear();
    use_delta_flag.clear();
    num_negative_pics = 0;
    num_positive_pics = 0;
    delta_poc_s0_minus1.clear();
    used_by_curr_pic_s0_flag.clear();
    delta_poc_s1_minus1.clear();
    used_by_curr_pic_s1_flag.clear();
}



void VuiParameters::toDefault()
{
    aspect_ratio_info_present_flag = 0;
    aspect_ratio_idc = 0;
    sar_width = 0;
    sar_height = 0;
    overscan_info_present_flag = 0;
    overscan_appropriate_flag = 0;
    video_signal_type_present_flag = 0;
    video_format = 5;
    video_full_range_flag = 0;
    colour_description_present_flag = 0;
    colour_primaries = 2;
    transfer_characteristics = 2;
    matrix_coeffs = 2;
    chroma_loc_info_present_flag = 0;
    chroma_sample_loc_type_top_field = 0;
    chroma_sample_loc_type_bottom_field = 0;
    neutral_chroma_indication_flag = 0;
    field_seq_flag = 0;
    frame_field_info_present_flag = 0;
    default_display_window_flag = 0;
    def_disp_win_left_offset = 0;
    def_disp_win_right_offset = 0;
    def_disp_win_top_offset = 0;
    def_disp_win_bottom_offset = 0;
    vui_timing_info_present_flag = 0;
    vui_num_units_in_tick = 0;
    vui_time_scale = 0;
    vui_poc_proportional_to_timing_flag = 0;
    vui_num_ticks_poc_diff_one_minus1 = 0;
    vui_hrd_parameters_present_flag = 0;
    hrd_parameters.toDefault();
    bitstream_restriction_flag = 0;
    tiles_fixed_structure_flag = 0;
    motion_vectors_over_pic_boundaries_flag = 0;
    restricted_ref_pic_lists_flag = 0;
    min_spatial_segmentation_idc = 0;
    max_bytes_per_pic_denom = 2;
    max_bits_per_min_cu_denom = 1;
    log2_max_mv_length_horizontal = 15;
    log2_max_mv_length_vertical = 15;
}

void VpsExtension::toDefault()
{
    splitting_flag = 0;
    scalability_mask_flag.clear();

    dimension_id_len_minus1.clear();
    vps_nuh_layer_id_present_flag = 0;
    layer_id_in_nuh.clear();
    dimension_id.clear();
    view_id_len = 0;
    view_id_val.clear();
    direct_dependency_flag.clear();
    num_add_layer_sets = 0;
    highest_layer_idx_plus1.clear();
    vps_sub_layers_max_minus1_present_flag = 0;
    sub_layers_vps_max_minus1.clear();
    max_tid_ref_present_flag = 0;
    max_tid_il_ref_pics_plus1.clear();
    default_ref_layers_active_flag = 0;
    vps_num_profile_tier_level_minus1 = 0;
    vps_profile_present_flag.clear();
    //profile_tier_level( vps_profile_present_flag[ i ], vps_max_sub_layers_minus1 )
    num_add_olss = 0;
    default_output_layer_idc = 0;
    //uint32_t NumOutputLayerSets() {  return num_add_olss + NumLayerSets();  }
    layer_set_idx_for_ols_minus1.clear();
    output_layer_flag.clear();
    profile_tier_level_idx.clear();
    alt_output_layer_flag.clear();
    vps_num_rep_formats_minus1 = 0;
    //rep_format
    rep_format_idx_present_flag = 0;
    vps_rep_format_idx.clear();
    max_one_active_ref_layer_flag = 0;
    vps_poc_lsb_aligned_flag = 0;
    poc_lsb_not_present_flag.clear();
    //dpb_size
    direct_dep_type_len_minus2 = 0;
    direct_dependency_all_layers_flag = 0;
    direct_dependency_all_layers_type = 0;
    direct_dependency_type.clear();
    vps_non_vui_extension_length = 0;
    vps_non_vui_extension_data_byte.clear();
    vps_vui_present_flag = 0;
    //alignbyte
    //vps_vui.toDefault();
}

void VpsVuiBspHrdParams::toDefault()
{
    vps_num_add_hrd_params = 0;
    cprms_add_present_flag.clear();
    num_sub_layer_hrd_minus1.clear();
    hrd_parameters.clear();
    num_signalled_partitioning_schemes.clear();
    num_partitions_in_scheme_minus1.clear();
    layer_included_in_partition_flag.clear();
    num_bsp_schedules_minus1.clear();
    bsp_hrd_idx.clear();
    bsp_sched_idx.clear();
}

void VpsVui::toDefault()
{
    cross_layer_pic_type_aligned_flag = 0;
    cross_layer_irap_aligned_flag = 0;
    all_layers_idr_aligned_flag = 0;
    bit_rate_present_vps_flag = 0;
    pic_rate_present_vps_flag = 0;
    bit_rate_present_flag.clear();
    pic_rate_present_flag.clear();
    avg_bit_rate.clear();
    max_bit_rate.clear();
    constant_pic_rate_idc.clear();
    avg_pic_rate.clear();

    video_signal_info_idx_present_flag =0;
    vps_num_video_signal_info_minus1 =0;

    video_signal_info.clear();
    vps_video_signal_info_idx.clear();
    tiles_not_in_use_flag = 0;
    tiles_in_use_flag.clear();
    loop_filter_not_across_tiles_flag.clear();
    tile_boundaries_aligned_flag.clear();

    wpp_not_in_use_flag = 0;
    wpp_in_use_flag.clear();
    single_layer_for_non_irap_flag = 0;
    higher_layer_irap_skip_flag = 0;
    ilp_restricted_ref_layers_flag = 0;
    min_spatial_segment_offset_plus1.clear();
    ctu_based_offset_enabled_flag.clear();
    min_horizontal_ctu_offset_plus1.clear();
    vps_vui_bsp_hrd_present_flag = 0;

    //vps_vui_bsp_hrd_params( )
    base_layer_parameter_set_compatibility_flag.clear();
}

void DpbSize::toDefault()
{
    sub_layer_flag_info_present_flag.clear();
    sub_layer_dpb_info_present_flag.clear();
    max_vps_dec_pic_buffering_minus1.clear();
    max_vps_num_reorder_pics.clear();
    max_vps_latency_increase_plus1.clear();
}

void RepFormat::toDefault()
{
    pic_width_vps_in_luma_samples = 0;
    pic_height_vps_in_luma_samples = 0;
    chroma_and_bit_depth_vps_present_flag =0;
    chroma_format_vps_idc = 0;
    separate_colour_plane_vps_flag = 0;
    bit_depth_vps_luma_minus8 = 0;
    bit_depth_vps_chroma_minus8 = 0;
    conformance_window_vps_flag = 0;
    conf_win_vps_left_offset = 0;
    conf_win_vps_right_offset = 0;
    conf_win_vps_top_offset = 0;
    conf_win_vps_bottom_offset = 0;
}

void VPS::toDefault()
{
  vps_video_parameter_set_id = 0;
  vps_max_layers_minus1 = 0;
  vps_max_sub_layers_minus1 = 0;
  vps_temporal_id_nesting_flag = 0;
  profile_tier_level.toDefault();
  vps_sub_layer_ordering_info_present_flag = 0;
  vps_max_dec_pic_buffering_minus1.clear();
  vps_max_num_reorder_pics.clear();
  vps_max_latency_increase_plus1.clear();
  vps_max_layer_id = 0;
  vps_num_layer_sets_minus1 = 0;
  layer_id_included_flag.clear();
  vps_timing_info_present_flag = 0;
  vps_num_units_in_tick = 0;
  vps_time_scale = 0;
  vps_poc_proportional_to_timing_flag = 0;
  vps_num_ticks_poc_diff_one_minus1 = 0;
  vps_num_hrd_parameters = 0;
  hrd_layer_set_idx.clear();
  cprms_present_flag.clear();
  hrd_parameters.clear();
  vps_extension_flag = 0;
  vps_extension.toDefault();

  vps_extension2_flag = 0;
  vps_3d_extension_flag = 0;
  vps_3d_extension.toDefault();

  vps_extension3_flag = 0;

}


void SPS::toDefault()
{
  sps_video_parameter_set_id = 0;
  sps_max_sub_layers_minus1 = 0;
  sps_temporal_id_nesting_flag = 0;
  profile_tier_level.toDefault();
  sps_seq_parameter_set_id = 0;
  chroma_format_idc = 0;
  separate_colour_plane_flag = 0;
  pic_width_in_luma_samples = 0;
  pic_height_in_luma_samples = 0;
  conformance_window_flag = 0;
  conf_win_left_offset = 0;
  conf_win_right_offset = 0;
  conf_win_top_offset = 0;
  conf_win_bottom_offset = 0;
  bit_depth_luma_minus8 = 0;
  bit_depth_chroma_minus8 = 0;
  log2_max_pic_order_cnt_lsb_minus4 = 0;
  sps_sub_layer_ordering_info_present_flag = 0;
  sps_max_dec_pic_buffering_minus1.clear();
  sps_max_num_reorder_pics.clear();
  sps_max_latency_increase_plus1.clear();
  log2_min_luma_coding_block_size_minus3 = 0;
  log2_diff_max_min_luma_coding_block_size = 0;
  log2_min_transform_block_size_minus2 = 0;
  log2_diff_max_min_transform_block_size = 0;
  max_transform_hierarchy_depth_inter = 0;
  max_transform_hierarchy_depth_intra = 0;
  scaling_list_enabled_flag = 0;
  scaling_list_data.toDefault();
  sps_scaling_list_data_present_flag = 0;
  amp_enabled_flag = 0;
  sample_adaptive_offset_enabled_flag = 0;
  pcm_enabled_flag = 0;
  pcm_sample_bit_depth_luma_minus1 = 0;
  pcm_sample_bit_depth_chroma_minus1 = 0;
  log2_min_pcm_luma_coding_block_size_minus3 = 0;
  log2_diff_max_min_pcm_luma_coding_block_size = 0;
  pcm_loop_filter_disabled_flag = 0;
  num_short_term_ref_pic_sets = 0;
  short_term_ref_pic_set.clear();
  long_term_ref_pics_present_flag = 0;
  num_long_term_ref_pics_sps = 0;
  lt_ref_pic_poc_lsb_sps.clear();
  used_by_curr_pic_lt_sps_flag.clear();
  sps_temporal_mvp_enabled_flag = 0;
  strong_intra_smoothing_enabled_flag = 0;
  vui_parameters_present_flag = 0;
  vui_parameters.toDefault();
  sps_extension_flag = 0;

  sps_range_extension_flag = 0;
  sps_multilayer_extension_flag = 0;
  sps_3d_extension_flag = 0;
  sps_extension_5bits = 0;

  sps_range_extension.toDefault();
  inter_view_mv_vert_constraint_flag = 0;
  sps_3d_extension.toDefault();
}



void PPS::toDefault()
{
    pps_pic_parameter_set_id = 0;
    pps_seq_parameter_set_id = 0;
    dependent_slice_segments_enabled_flag = 0;
    output_flag_present_flag = 0;
    num_extra_slice_header_bits = 0;
    sign_data_hiding_flag = 0;
    cabac_init_present_flag = 0;
    num_ref_idx_l0_default_active_minus1 = 0;
    num_ref_idx_l1_default_active_minus1 = 0;
    init_qp_minus26 = 0;
    constrained_intra_pred_flag = 0;
    transform_skip_enabled_flag = 0;
    cu_qp_delta_enabled_flag = 0;
    diff_cu_qp_delta_depth = 0;
    pps_cb_qp_offset = 0;
    pps_cr_qp_offset = 0;
    pps_slice_chroma_qp_offsets_present_flag = 0;
    weighted_pred_flag = 0;
    weighted_bipred_flag = 0;
    transquant_bypass_enabled_flag = 0;
    tiles_enabled_flag = 0;
    entropy_coding_sync_enabled_flag = 0;
    num_tile_columns_minus1 = 0;
    num_tile_rows_minus1 = 0;
    uniform_spacing_flag = 1;
    column_width_minus1.clear();
    row_height_minus1.clear();
    loop_filter_across_tiles_enabled_flag = 0;
    pps_loop_filter_across_slices_enabled_flag = 0;
    deblocking_filter_control_present_flag = 0;
    deblocking_filter_override_enabled_flag = 0;
    pps_deblocking_filter_disabled_flag = 0;
    pps_beta_offset_div2 = 0;
    pps_tc_offset_div2 = 0;
    pps_scaling_list_data_present_flag = 0;
    lists_modification_present_flag = 0;
    log2_parallel_merge_level_minus2 = 0;
    slice_segment_header_extension_present_flag = 0;
    pps_extension_flag = 0;

    pps_range_extension_flag = 0;
    pps_multilayer_extension_flag = 0;
    pps_3d_extension_flag = 0;
    pps_extension_5bits = 0;

    pps_range_extension.toDefault();

    pps_multilayer_extension.toDefault();

    pps_3d_extension.toDefault();
}



void Slice::toDefault()
{
  first_slice_segment_in_pic_flag = 0;
  no_output_of_prior_pics_flag = 0;
  slice_pic_parameter_set_id = 0;
  dependent_slice_segment_flag = 0;
  slice_segment_address = 0;
  slice_reserved_undetermined_flag.clear();
  slice_type = 0;
  pic_output_flag = 1;
  colour_plane_id = 0;
  slice_pic_order_cnt_lsb = 0;
  short_term_ref_pic_set_sps_flag = 0;
  short_term_ref_pic_set_idx = 0;
  num_ref_idx_l0_active_minus1 = 0;
  num_ref_idx_l1_active_minus1 = 0;
  slice_temporal_mvp_enabled_flag = 0;
  collocated_from_l0_flag = 1;
  deblocking_filter_override_flag = 0;
  slice_sao_luma_flag = 0;
  slice_sao_chroma_flag = 0;
  num_long_term_pics = 0;
  num_long_term_sps = 0;
}


void ScalingListData::toDefault()
{
  scaling_list_delta_coef.clear();
  scaling_list_pred_mode_flag.clear();
  scaling_list_pred_matrix_id_delta.clear();
  scaling_list_dc_coef_minus8.clear();
}


void AUD::toDefault()
{
  pic_type = 0;
}


void SEI::toDefault()
{
  sei_message.clear();
}


void SeiMessage::toDefault()
{
  num_payload_type_ff_bytes = 0;
  num_payload_size_ff_bytes = 0;
  last_payload_type_byte = 0;
  last_payload_size_byte = 0;
}


void DecodedPictureHash::toDefault()
{
  hash_type = 0;
}

void RefPicListModification::toDefault()
{
  ref_pic_list_modification_flag_l0 = 0;
  list_entry_l0.clear();
  ref_pic_list_modification_flag_l1 = 0;
  list_entry_l1.clear();
}



void PredWeightTable::toDefault()
{
  luma_log2_weight_denom = 0;
  delta_chroma_log2_weight_denom = 0;

  luma_weight_l0_flag.clear();
  chroma_weight_l0_flag.clear();
  delta_luma_weight_l0.clear();
  luma_offset_l0.clear();
  delta_chroma_weight_l0.clear();
  delta_chroma_offset_l0.clear();
  luma_weight_l1_flag.clear();
  chroma_weight_l1_flag.clear();
  delta_luma_weight_l1.clear();
  luma_offset_l1.clear();
  delta_chroma_weight_l1.clear();
  delta_chroma_offset_l1.clear();
}


void MasteringDisplayInfo::toDefault()
{
  memset(display_primary_x, 0, sizeof(display_primary_x));
  memset(display_primary_y, 0, sizeof(display_primary_x));
  white_point_x = 0;
  white_point_y = 0;
  max_display_mastering_luminance = 0;
  min_display_mastering_luminance = 0;
}


void BufferingPeriod::toDefault()
{
  nal_initial_cpb_removal_delay.clear();
  nal_initial_cpb_removal_offset.clear();
  nal_initial_alt_cpb_removal_delay.clear();
  nal_initial_alt_cpb_removal_offset.clear();
  vcl_initial_cpb_removal_delay.clear();
  vcl_initial_cpb_removal_offset.clear();
  vcl_initial_alt_cpb_removal_delay.clear();
  vcl_initial_alt_cpb_removal_offset.clear();

  irap_cpb_params_present_flag = 0;
  cpb_delay_offset = 0;
  dpb_delay_offset = 0;

}


void PicTiming::toDefault()
{
  pic_struct = 0;
  source_scan_type = 0;
  duplicate_flag = 0;
  au_cpb_removal_delay_minus1 = 0;
  pic_dpb_output_delay = 0;
  pic_dpb_output_du_delay = 0;
  num_decoding_units_minus1 = 0;
  du_common_cpb_removal_delay_flag = 0;
  du_common_cpb_removal_delay_increment_minus1 = 0;
  num_nalus_in_du_minus1.clear();
  du_cpb_removal_delay_increment_minus1.clear();
}
