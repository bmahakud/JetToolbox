###############################################
####
####   Jet Substructure Toolbox (jetToolBox)
####   Python function for easy access of 
####   jet substructure tools implemented in CMS
####
####   Alejandro Gomez Espinosa (gomez@physics.rutgers.edu)
####
###############################################
import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSSoftDropMass, ak8PFJetsCHSPruned, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSTrimmed, ak8PFJetsCHSTrimmedMass, ak8PFJetsCHSFiltered, ak8PFJetsCHSFilteredMass, ak4PFJetsCHS, ca15PFJetsCHSMassDropFiltered, hepTopTagPFJetsCHS, ak8PFJetsCHSConstituents, puppi
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection, updateJetCollection


def jetToolbox( proc, jetType, jetSequence, outputFile,
		updateCollection='', updateCollectionSubjets='',
		newPFCollection=False, nameNewPFCollection = '',	
		PUMethod='CHS', 
		miniAOD=True,
		runOnMC=True,
		JETCorrPayload='', JETCorrLevels = [ 'None' ], GetJetMCFlavour=True,
		Cut = '', 
		postFix='',
		bTagDiscriminators = None, 
		bTagInfos = None, 
		subJETCorrPayload='', subJETCorrLevels = [ 'None' ], GetSubjetMCFlavour=True,
		CutSubjet = '', 
		addPruning=False, zCut=0.1, rCut=0.5, addPrunedSubjets=False,
		addSoftDrop=False, betaCut=0.0,  zCutSD=0.1, addSoftDropSubjets=False,
		addTrimming=False, rFiltTrim=0.2, ptFrac=0.03,
		addFiltering=False, rfilt=0.3, nfilt=3,
		addCMSTopTagger=False,
		addMassDrop=False,
		addHEPTopTagger=False,
		addNsub=False, maxTau=4,
		addNsubSubjets=False, subjetMaxTau=4,
		addPUJetID=False,
		addQJets=False,
		addQGTagger=False, QGjetsLabel='chs',
		addEnergyCorrFunc=False, ecfBeta = 1.0, maxECF=3,
		):

	runOnData = not runOnMC
	if runOnData:
		GetJetMCFlavour = False
		GetSubjetMCFlavour = False
	
	###############################################################################
	#######  Verifying some inputs and defining variables
	###############################################################################
	print '|---- jetToolbox: Initialyzing collection...'
	if newPFCollection: print '|---- jetToolBox: Using '+ nameNewPFCollection +' as PFCandidates collection'
	supportedJetAlgos = { 'ak': 'AntiKt', 'ca' : 'CambridgeAachen', 'kt' : 'Kt' }
	recommendedJetAlgos = [ 'ak4', 'ak8', 'ca4', 'ca8', 'ca10', 'ca15' ]
	payloadList = [ 'None',
			'AK1PFchs', 'AK2PFchs', 'AK3PFchs', 'AK4PFchs', 'AK5PFchs', 'AK6PFchs', 'AK7PFchs', 'AK8PFchs', 'AK9PFchs', 'AK10PFchs',
			'AK1PFPuppi', 'AK2PFPuppi', 'AK3PFPuppi', 'AK4PFPuppi', 'AK5PFPuppi', 'AK6PFPuppi', 'AK7PFPuppi', 'AK8PFPuppi', 'AK9PFPuppi', 'AK10PFPuppi',  
			'AK1PFSK', 'AK2PFSK', 'AK3PFSK', 'AK4PFSK', 'AK5PFSK', 'AK6PFSK', 'AK7PFSK', 'AK8PFSK', 'AK9PFSK', 'AK10PFSK',  
			'AK1PF', 'AK2PF', 'AK3PF', 'AK4PF', 'AK5PF', 'AK6PF', 'AK7PF', 'AK8PF', 'AK9PF', 'AK10PF' ]
	JECLevels = [ 'L1Offset', 'L1FastJet', 'L1JPTOffset', 'L2Relative', 'L3Absolute', 'L5Falvour', 'L7Parton' ]
	if runOnData:
		JECLevels += ['L2L3Residual']
	jetAlgo = ''
	algorithm = ''
	size = ''
	for TYPE, tmpAlgo in supportedJetAlgos.iteritems(): 
		if TYPE in jetType.lower():
			jetAlgo = TYPE
			algorithm = tmpAlgo
			size = jetType.replace( TYPE, '' )

	jetSize = 0.
	if int(size) in range(0, 20): jetSize = int(size)/10.
	else: print '|---- jetToolBox: jetSize has not a valid value. Insert a number between 1 and 20 after algorithm, like: AK8'
	### Trick for uppercase/lowercase algo name
        
	jetALGO = jetAlgo.upper()+size
	jetalgo = jetAlgo.lower()+size
        print 'jetAlgo: ',jetAlgo
        print 'jetALGO: ',jetALGO
        print 'jetalgo', jetalgo
	if jetalgo not in recommendedJetAlgos: print '|---- jetToolBox: CMS recommends the following jet algoritms: '+' '.join(recommendedJetAlgos)+'. You are using', jetalgo,'.'


	#################################################################################
	####### Toolbox start 
	#################################################################################

	elemToKeep = []
	jetSeq = cms.Sequence()
	genParticlesLabel = ''
	pvLabel = ''
	tvLabel = ''
	toolsUsed = []

	### List of Jet Corrections
	if not set(JETCorrLevels).issubset(set(JECLevels)): 
		if ( 'CHS' in PUMethod ) or  ( 'Plain' in PUMethod ): JETCorrLevels = ['L1FastJet','L2Relative', 'L3Absolute']
		else: JETCorrLevels = [ 'L2Relative', 'L3Absolute']
		if runOnData: JETCorrLevels.append('L2L3Residual')
	if not set(subJETCorrLevels).issubset(set(JECLevels)): 
		if ( 'CHS' in PUMethod ) or  ( 'Plain' in PUMethod ): subJETCorrLevels = ['L1FastJet','L2Relative', 'L3Absolute']
		else: subJETCorrLevels = [ 'L2Relative', 'L3Absolute']
		if runOnData: subJETCorrLevels.append('L2L3Residual')


	if not updateCollection: 

		## b-tag discriminators
		if bTagDiscriminators is None:
			bTagDiscriminators = [
					'pfTrackCountingHighEffBJetTags',
					'pfTrackCountingHighPurBJetTags',
					'pfJetProbabilityBJetTags',
					'pfJetBProbabilityBJetTags',
					'pfSimpleSecondaryVertexHighEffBJetTags',
					'pfSimpleSecondaryVertexHighPurBJetTags',
					'pfCombinedSecondaryVertexV2BJetTags',
					'pfCombinedInclusiveSecondaryVertexV2BJetTags',
					'pfCombinedMVAV2BJetTags'
			]

		#### For MiniAOD
		if miniAOD:

			print '|---- jetToolBox: JETTOOLBOX RUNNING ON MiniAOD FOR '+jetALGO+' JETS USING '+PUMethod

			genParticlesLabel = 'prunedGenParticles'
			pvLabel = 'offlineSlimmedPrimaryVertices'
			svLabel = 'slimmedSecondaryVertices'
			tvLabel = 'unpackedTracksAndVertices'
			muLabel = 'slimmedMuons'
			elLabel = 'slimmedElectrons'
			pfCand =  'packedPFCandidates'

			if runOnMC:
				## Filter out neutrinos from packed GenParticles
				setattr( proc, 'packedGenParticlesForJetsNoNu', 
						cms.EDFilter("CandPtrSelector", 
							src = cms.InputTag("packedGenParticles"), 
							cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
							))
				jetSeq += getattr(proc, 'packedGenParticlesForJetsNoNu' )

				setattr( proc, jetalgo+'GenJetsNoNu', 
						ak4GenJets.clone( src = 'packedGenParticlesForJetsNoNu', 
							rParam = jetSize, 
							jetAlgorithm = algorithm ) ) 
				jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )

			#for Inclusive Vertex Finder
			proc.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

			

		####  Creating PATjets
		tmpPfCandName = pfCand.lower()
		if 'Puppi' in PUMethod:
			if ('puppi' in tmpPfCandName): 
				srcForPFJets = pfCand
				print '|---- jetToolBox: Not running puppi algorithm because keyword puppi was specified in nameNewPFCollection, but applying puppi corrections.' 
			else: 
				proc.load('CommonTools.PileupAlgos.Puppi_cff')
				puppi.candName = cms.InputTag( pfCand ) 
				if miniAOD:
				  puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
				  puppi.clonePackedCands = cms.bool(True)
				jetSeq += getattr(proc, 'puppi' )
				srcForPFJets = 'puppi'

			from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
			setattr( proc, jetalgo+'PFJetsPuppi', 
					ak4PFJetsPuppi.clone( src = cms.InputTag( srcForPFJets ),
						doAreaFastjet = True, 
						rParam = jetSize, 
						jetAlgorithm = algorithm ) )  
			jetSeq += getattr(proc, jetalgo+'PFJetsPuppi' )

			#if JETCorrPayload not in payloadList: JETCorrPayload = 'AK'+size+'PFPuppi'
			#if subJETCorrPayload not in payloadList: subJETCorrPayload = 'AK4PFPuppi'
		if 'None' in JETCorrPayload: JEC = None
		#else: JEC = ( JETCorrPayload.replace('CS','chs').replace('SK','chs') , JETCorrLevels, 'None' )   ### temporary
		#else: JEC = ( JETCorrPayload., JETCorrLevels, 'None' ) 
		#print '|---- jetToolBox: Applying this corrections: '+str(JEC)

		if ( PUMethod in [ 'CHS', 'CS', 'Puppi' ] ) and miniAOD: setattr( proc, jetalgo+'PFJets'+PUMethod+'Constituents', cms.EDProducer("MiniAODJetConstituentSelector", src = cms.InputTag( jetalgo+'PFJets'+PUMethod ), cut = cms.string( Cut ) ))
	#	else: setattr( proc, jetalgo+'PFJets'+PUMethod+'Constituents', cms.EDProducer("PFJetConstituentSelector", src = cms.InputTag( jetalgo+'PFJets'+PUMethod ), cut = cms.string( Cut ) ))
	#	jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'Constituents' )
                
		addJetCollection(
				proc,
				labelName = jetALGO+'PF'+PUMethod,
				jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod),
				postfix = postFix, 
				algo = jetalgo,
				rParam = jetSize,
				jetCorrections = JEC if JEC is not None else None, 
				pfCandidates = cms.InputTag( pfCand ),  
				svSource = cms.InputTag( svLabel ),  
				genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
				pvSource = cms.InputTag( pvLabel ), 
				muSource = cms.InputTag( muLabel ),
				elSource = cms.InputTag( elLabel ),
				btagDiscriminators = bTagDiscriminators,
				btagInfos = bTagInfos,
				getJetMCFlavour = GetJetMCFlavour,
				genParticles = cms.InputTag(genParticlesLabel),
				outputModules = ['outputFile']
				) 
		
		#getattr(proc,'patJets'+jetALGO+'PF'+PUMethod+postFix).addTagInfos = cms.bool(True)

		#if 'CS' in PUMethod: getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod+postFix ).getJetMCFlavour = False  # CS jets cannot be re-clustered from their constituents
		#patJets = 'patJets'
		#patSubJets = ''
		
             
