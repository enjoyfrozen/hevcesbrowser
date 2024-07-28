#ifndef HEVC_PARSER_IMPL
#define HEVC_PARSER_IMPL



#include "HevcParser.h"
#include "BitstreamReader.h"

#include <map>
#include <list>
#include <memory>

namespace HEVC
{
  class HevcParserImpl: public Parser
  {
    public:
      virtual std::size_t process(const uint8_t *pdata, std::size_t size, std::size_t offset = 0);

      virtual void addConsumer(Consumer *pconsumer);
      virtual void releaseConsumer(Consumer *pconsumer);

    protected:
      void processNALUnit(const uint8_t *pdata, std::size_t size, const Parser::Info &info);
      void processNALUnitHeader(BitstreamReader &bs, NALHeader* header);
      void processVPS(std::shared_ptr<VPS> pvps, BitstreamReader &bs, const Parser::Info &info);
      void processSPS(std::shared_ptr<SPS> psps, BitstreamReader &bs, const Parser::Info &info);
      void processPPS(std::shared_ptr<PPS> ppps, BitstreamReader &bs, const Parser::Info &info);
      void processSlice(std::shared_ptr<Slice> pslice, BitstreamReader &bs, const Parser::Info &info);
      void processAUD(std::shared_ptr<AUD> paud, BitstreamReader &bs, const Parser::Info &info);


      void processSEI(std::shared_ptr<SEI> psei, BitstreamReader &bs, const Parser::Info &info);
      void processSliceHeader(std::shared_ptr<Slice> pslice, BitstreamReader &bs, const Parser::Info &info);
      void processSliceData(std::shared_ptr<Slice> pslice, BitstreamReader &bs, const Parser::Info &info);
      ProfileTierLevel processProfileTierLevel(uint8_t profilePresentFlag, std::size_t max_sub_layers_minus1, BitstreamReader &bs, const Parser::Info &info);
      HrdParameters processHrdParameters(uint8_t commonInfPresentFlag, std::size_t maxNumSubLayersMinus1, BitstreamReader &bs);
      SubLayerHrdParameters processSubLayerHrdParameters(uint8_t sub_pic_hrd_params_present_flag, std::size_t CpbCnt, BitstreamReader &bs);
      ShortTermRefPicSet processShortTermRefPicSet(std::size_t stRpsIdx, size_t num_short_term_ref_pic_sets, const std::vector<ShortTermRefPicSet> &refPicSets, std::shared_ptr<SPS> psps, BitstreamReader &bs, const Parser::Info &info);
      VuiParameters processVuiParameters(std::size_t sps_max_sub_layers_minus1, BitstreamReader &bs);
      ScalingListData processScalingListData(BitstreamReader &bs);
      RefPicListModification processRefPicListModification(BitstreamReader &bs, std::shared_ptr<Slice> pslice);
      PredWeightTable processPredWeightTable(BitstreamReader &bs, std::shared_ptr<Slice> pslice);
      VpsVui processVpsVui(std::shared_ptr<VPS> pvps, BitstreamReader& bs);
      VpsVuiBspHrdParams processVpsVuiBspHrdParams(std::shared_ptr<VPS> pvps, BitstreamReader& bs);

      //zorro add
      int NumOutputLayerSets = 0;
      std::vector<uint8_t> NumLayersInIdList;
      //(F-10£©
      std::vector<uint32_t> MaxSubLayersInLayerSetMinus1;
      void updateMaxSubLayersInLayerSetMinus1(std::shared_ptr<VPS> pvps);
      
      //(F-3)
      int NumViews = 0;
      std::vector < std::vector<uint8_t>> ScalabilityId;
      std::vector<uint32_t>  DepthLayerFlag, ViewOrderIdx, DependencyId, AuxId;
      std::vector<uint8_t> LayerIdxInVps;
      void processNumViews(std::shared_ptr<VPS> pvps);

      //(F-4)
      std::vector < std::vector<uint8_t>> DependencyFlag;
      void updateDependencyFlag(std::shared_ptr<VPS> pvps);

      //(F-5)
      std::vector < std::vector<uint8_t>> IdDirectRefLayer, IdRefLayer, IdPredictedLayer;
      std::vector<uint8_t> NumDirectRefLayers, NumRefLayers, NumPredictedLayers;
      void updateF5(std::shared_ptr<VPS> pvps);

      //(F-6) 
      int NumIndependentLayers = 0;
      std::vector<std::vector<uint8_t>> TreePartitionLayerIdList;
      std::vector<uint8_t> NumLayersInTreePartition;
      void updateNumIndependentLayers(std::shared_ptr<VPS> pvps);


      //(7-3)
      int MaxLayersMinus1 = 0;
      int NumScalabilityTypes = 0;
      void updateMaxLayersMinus1(std::shared_ptr<VPS> pvps);

      //(F-7)
      int NumLayerSets = 0;
      void updateNumLayerSets(std::shared_ptr<VPS> pvps);

       //(F-8)
      int FirstAddLayerSetIdx = 0;
      int LastAddLayerSetIdx = 0;
      void updateF8(std::shared_ptr<VPS> pvps);

      //(F-9) update LayerSetLayerIdList, and NumLayersInIdList
      std::vector<std::vector<uint16_t>> LayerSetLayerIdList;

      std::vector<std::vector<uint16_t>> OutputLayerFlag;
      void updateNumLayersInIdList(int i, std::shared_ptr<VPS> pvps);

      // The nuh_layer_id value of the NAL unit containing the PPS that is activated for a layer layerA with nuh_layer_id equal to nuhLayerIdA shall be equal to 0, or nuhLayerIdA, or the nuh_layer_id of a direct or indirect reference layer of layerA.
      int vps_nuhLayerIdA = 0; //zorro: default 0
      int sps_nuhLayerIdA = 0; //zorro: default 0
      int pps_nuhLayerIdA = 0; //zorro: default 0
      void updateOutputLayerFlag(std::shared_ptr<VPS> pvps);

      //(F-11)
      std::vector<uint8_t> OlsIdxToLsIdx;
      void ProcessOlsIdxToLsIdx(int maxIdx, std::shared_ptr<VPS> pvps);
      
      //(F-12)
      int defaultOutputLayerIdc = 0;
      std::vector<uint8_t> NumOutputLayersInOutputLayerSet, OlsHighestOutputLayerId;
      void ProcessNumOutputLayersInOutputLayerSet(std::shared_ptr<VPS> pvps);
      
      //(F-13)
      std::vector<std::vector<uint8_t>> NecessaryLayerFlag;
      std::vector<uint8_t> NumNecessaryLayers;
      void ProcessF13();

      DpbSize processDpbSize(std::shared_ptr<VPS> pvps, BitstreamReader& bs);
      RepFormat processRepFormat(BitstreamReader& bs);

      VideoSignalInfo processVideoSignalInfo(BitstreamReader& bs);

      void ProcessVpsExtension(std::shared_ptr<VPS> pvps, BitstreamReader& bs);

      //add for sps
      SpsRangeExtension processSpsRangeExtension(BitstreamReader& bs);
      Sps3dExtension processSps3dExtension(BitstreamReader& bs);
      //add end

      //add for pps
      PpsRangeExtension processPpsRangeExtension(std::shared_ptr<PPS> ppps, BitstreamReader &bs);
      PpsMultilayerExtension processPpsMultilayerExtension(std::shared_ptr<PPS> ppps, BitstreamReader& bs);
      Pps3dExtension processPps3dExtension(std::shared_ptr<PPS> ppps, BitstreamReader& bs);
      DeltaDlt processDeltaDlt(int pps_bit_depth_for_depth_layers_minus8, BitstreamReader& bs);
      ColourMappingTable processColourMappingTable(BitstreamReader& bs);
      //add end

      void processDecodedPictureHash(std::shared_ptr<DecodedPictureHash> pdecPicHash, BitstreamReader &bs);
      void processUserDataUnregistered(std::shared_ptr<UserDataUnregistered> pSeiPayload, BitstreamReader &bs, std::size_t payloadSize);
      void processSceneInfo(std::shared_ptr<SceneInfo> pSeiPayload, BitstreamReader &bs);
      void processFullFrameSnapshot(std::shared_ptr<FullFrameSnapshot> pSeiPayload, BitstreamReader &bs);
      void processProgressiveRefinementSegmentStart(std::shared_ptr<ProgressiveRefinementSegmentStart> pSeiPayload, BitstreamReader &bs);
      void processProgressiveRefinementSegmentEnd(std::shared_ptr<ProgressiveRefinementSegmentEnd> pSeiPayload, BitstreamReader &bs);
      void processBufferingPeriod(std::shared_ptr<BufferingPeriod> pSeiPayload, BitstreamReader &bs);
      void processPicTiming(std::shared_ptr<PicTiming> pSeiPayload, BitstreamReader &bs);
      void processRecoveryPoint(std::shared_ptr<RecoveryPoint> pSeiPayload, BitstreamReader &bs);
      void processActiveParameterSets(std::shared_ptr<ActiveParameterSets> pSeiPayload, BitstreamReader &bs);
      void processTemporalLevel0Index(std::shared_ptr<TemporalLevel0Index> pSeiPayload, BitstreamReader &bs);
      void processRegionRefreshInfo(std::shared_ptr<RegionRefreshInfo> pSeiPayload, BitstreamReader &bs);
      void processToneMapping(std::shared_ptr<ToneMapping> pSeiPayload, BitstreamReader &bs);
      void processFramePacking(std::shared_ptr<FramePacking> pSeiPayload, BitstreamReader &bs);
      void processDisplayOrientation(std::shared_ptr<DisplayOrientation> pSeiPayload, BitstreamReader &bs);
      void processSOPDescription(std::shared_ptr<SOPDescription> pSeiPayload, BitstreamReader &bs);
      void processTimeCode(std::shared_ptr<TimeCode> pSeiPayload, BitstreamReader &bs);
      void processMasteringDisplayInfo(std::shared_ptr<MasteringDisplayInfo> pcolourVolume, BitstreamReader &bs);
      void processSegmRectFramePacking(std::shared_ptr<SegmRectFramePacking> pSeiPayload, BitstreamReader &bs);
      void processKneeFunctionInfo(std::shared_ptr<KneeFunctionInfo> pSeiPayload, BitstreamReader &bs);
      void processChromaResamplingFilterHint(std::shared_ptr<ChromaResamplingFilterHint> pSeiPayload, BitstreamReader &bs);
      void processColourRemappingInfo(std::shared_ptr<ColourRemappingInfo> pSeiPayload, BitstreamReader &bs);
      void processContentLightLevelInfo(std::shared_ptr<ContentLightLevelInfo> pSeiPayload, BitstreamReader &bs);
      void processAlternativeTransferCharacteristics(std::shared_ptr<AlternativeTransferCharacteristics> pSeiPayload, BitstreamReader &bs);

      void onWarning(const std::string &warning, const Info *pInfo, WarningType type);

      std::map<uint32_t, std::shared_ptr<VPS> >          m_vpsMap;
      std::map<uint32_t, std::shared_ptr<SPS> >          m_spsMap;
      std::map<uint32_t, std::shared_ptr<PPS> >          m_ppsMap;
      std::shared_ptr<Slice>                             m_lastSlice;

      std::list<Consumer *>          m_consumers;
  };
}

#endif
