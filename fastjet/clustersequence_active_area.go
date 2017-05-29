// Copyright 2017 The go-hep Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fastjet

import (
	"math"

	"go-hep.org/x/hep/fmom"
)

type ClusterSequenceActiveArea struct {
	cs        *ClusterSequence
	def       JetDefinition
	ghostSpec *GhostedAreaSpec

	ghostSpecRepeat int
	maxRapForArea   float64
	safeRapForArea  float64
	averageArea     []float64
	averageArea2    []float64
	averageP4       []Jet // FIXME averageP4 might not be needed
	nonJetArea      float64
	nonJetArea2     float64
	nonJetN         float64
}

func NewClusterSequenceActiveArea(jets []Jet, jetDefIn JetDefinition, ghostSpec *GhostedAreaSpec) (*ClusterSequenceActiveArea, error) {
	aa := &ClusterSequenceActiveArea{
		def:       jetDefIn,
		ghostSpec: ghostSpec,
	}
	var err error

	// FIXME check if ClusterSequence needs to be run here
	aa.cs, err = constructClusterSequence(jets, jetDefIn)
	if err != nil {
		return nil, err
	}
	return aa, nil
}

func (aa *ClusterSequenceActiveArea) initializeAndRunAA() error {
	continueRun, err := aa.initialize()
	if err != nil {
		return err
	}
	if continueRun {
		err := aa.run()
		if err != nil {
			return err
		}
		aa.postProcess()
	}
	return nil
}

func (aa *ClusterSequenceActiveArea) initialize() (bool, error) {
	aa.ghostSpecRepeat = aa.ghostSpec.repeat

	aa.averageArea = make([]float64, 2*len(aa.cs.jets))
	aa.averageArea2 = make([]float64, 2*len(aa.cs.jets))

	aa.resizeAndZero()
	aa.maxRapForArea = aa.ghostSpec.ghostMaxRap
	aa.safeRapForArea = aa.maxRapForArea - aa.def.R()

	if aa.ghostSpec.repeat <= 0 {
		err := aa.initializeAndRunAA()
		if err != nil {
			return false, err
		}
		return false, nil
	}

	// decant_options stores variables

	err := aa.cs.init()
	if err != nil {
		return false, err
	}

	return true, nil
}

func (aa *ClusterSequenceActiveArea) resizeAndZero() {
	for i := range aa.averageArea {
		aa.averageArea[i] = 0
	}
	for i := range aa.averageArea2 {
		aa.averageArea2[i] = 0
	}
	for i := range aa.averageP4 {
		aa.averageP4[i] = NewJet(0, 0, 0, 0)
	}
	aa.nonJetArea = 0
	aa.nonJetArea2 = 0
	aa.nonJetN = 0
}

func (aa *ClusterSequenceActiveArea) run() error {
	inputJets := aa.cs.jets

	uniqueTree := make([]int, 0)

	for irepeat := 0; irepeat < aa.ghostSpec.repeat; irepeat++ {
		clustSeq, err := NewClusterSequenceActiveAreaExplicitGhosts(inputJets, aa.def, aa.ghostSpec)
		if err != nil {
			return err
		}
		if irepeat == 0 {
			// returning clustSeq shouldn't be necessary
			var err error
			clustSeq, err = aa.transferGhostFreeHistory(clustSeq)
			if err != nil {
				return err
			}
			uniqueTree = aa.cs.uniqueHistoryOrder()
		}
		uniqueTree, clustSeq = aa.transferAreas(uniqueTree, clustSeq)
	}
	return nil
}

// return values should be unnecessary
func (aa *ClusterSequenceActiveArea) transferAreas(uniqueHist []int, gs *ClusterSequenceActiveAreaExplicitGhosts) ([]int, *ClusterSequenceActiveAreaExplicitGhosts) {
	gsHistory := gs.cs.history
	gsJets := gs.cs.jets
	gsUniqueHist := gs.cs.uniqueHistoryOrder()

	j := -1
	histInd := -1

	ourAreas := make([]float64, len(aa.cs.history))
	for i := range ourAreas {
		ourAreas[i] = 0
	}

	our4Vectors := make([]Jet, len(aa.cs.history))
	for i := range our4Vectors {
		our4Vectors[i] = NewJet(0, 0, 0, 0)
	}

	for i := 0; i < len(gsHistory); i++ {
		gsHistInd := gsUniqueHist[i]
		if gsHistInd < gs.cs.initn {
			continue
		}
		gsHist := gsHistory[gsUniqueHist[i]]
		parent1 := gsHist.parent1
		parent2 := gsHist.parent2

		if parent2 == beamJetIndex {
			jet := gsJets[gsHistory[parent1].jet]
			areaLocal := gs.areas[jet.hidx]
			extArea := gs.area4vectors[jet.hidx]

			if gs.isPureGhost[parent1] {
				// record ghost jets. Not sure if necessary
				// aa.ghosts = append(aa.ghosts,*jet.GhostJet(areaLocal))
				if math.Abs(jet.rap) < aa.safeRapForArea {
					aa.nonJetArea += areaLocal
					aa.nonJetArea2 += areaLocal * areaLocal
					aa.nonJetN += 1
				}
			} else {
				j++
				for ; j < len(aa.cs.history); j++ {
					histInd = uniqueHist[j]
					if histInd >= aa.cs.initn {
						break
					}
				}

				// left out check

				// refJetInd := aa.cs.history[aa.cs.history[histInd].parent1].jet

				// check

				// refJet := aa.cs.jets[refJetInd]

				// more checks

				ourAreas[histInd] = areaLocal
				our4Vectors[histInd] = extArea

				if aa.cs.history[histInd].parent1 != inexistentParent {
					ourAreas[aa.cs.history[histInd].parent1] = areaLocal
					our4Vectors[aa.cs.history[histInd].parent1] = extArea
				}
			}
		} else if !gs.isPureGhost[parent1] && !gs.isPureGhost[parent2] {
			j++
			for ; j < len(aa.cs.history); j++ {
				histInd = uniqueHist[j]
				if histInd >= aa.cs.initn {
					break
				}
			}

			jet := gsJets[gsHist.jet]
			// jetRef := aa.cs.jets[aa.cs.history[histInd].jet]

			// fmt.Println(jet.E()," ",jet.P2())
			// fmt.Println(jetRef.E()," ",jetRef.P2())
			// check

			areasLocal := gs.areas[jet.hidx]
			ourAreas[histInd] += areasLocal

			extArea := gs.area4vectors[jet.hidx]

			// what is this recombiner doing?
			var err error
			our4Vectors[histInd], err = aa.def.recombiner.Recombine(&our4Vectors[histInd], &extArea)
			if err != nil {

			}

			jet1 := gsJets[gsHistory[parent1].jet]
			ourParent1 := aa.cs.history[histInd].parent1
			if ourParent1 != inexistentParent {
				ourAreas[ourParent1] = gs.areas[jet1.hidx]
				our4Vectors[ourParent1] = gs.area4vectors[jet1.hidx]
			}
			jet2 := gsJets[gsHistory[parent2].jet]
			ourParent2 := aa.cs.history[histInd].parent2
			if ourParent2 != inexistentParent {
				ourAreas[ourParent2] = gs.areas[jet2.hidx]
				our4Vectors[ourParent2] = gs.area4vectors[jet2.hidx]
			}
		}
	}

	// missing stuff

	for areaInd := 0; areaInd < len(ourAreas); areaInd++ {
		aa.averageArea[areaInd] += ourAreas[areaInd]
		aa.averageArea2[areaInd] += ourAreas[areaInd] * ourAreas[areaInd]
	}

	for i := 0; i < len(our4Vectors); i++ {
		// recombine plus equal
	}

	return uniqueHist, gs
}

func (aa *ClusterSequenceActiveArea) transferGhostFreeHistory(gs *ClusterSequenceActiveAreaExplicitGhosts) (*ClusterSequenceActiveAreaExplicitGhosts, error) {
	gsHistory := gs.cs.history
	gsToSelfHistMap := make([]int, len(gsHistory))

	// transfer strategy used

	igs := 0
	iself := 0
	for ; igs < len(gsHistory) && gsHistory[igs].parent1 == inexistentParent; igs++ {
		if !gs.isPureGhost[igs] {
			gsToSelfHistMap[igs] = iself
			iself++
		} else {
			gsToSelfHistMap[igs] = invalidIndex
		}
	}

	if iself != len(aa.cs.history) {
		panic("iself doesn't match aa history")
	}

	if igs == len(gsHistory) {
		return gs, nil
	}

	for ; igs < len(gsHistory); igs++ {
		if gs.isPureGhost[igs] {
			gsToSelfHistMap[igs] = invalidIndex
			continue
		}
		histElem := gsHistory[igs]

		parent1isGhost := gs.isPureGhost[histElem.parent1]
		parent2isGhost := false
		if histElem.parent2 != beamJetIndex {
			parent2isGhost = gs.isPureGhost[histElem.parent2]
		}

		if parent1isGhost && !parent2isGhost && histElem.parent2 >= 0 {
			gsToSelfHistMap[igs] = gsToSelfHistMap[histElem.parent2]
			continue
		}
		if !parent1isGhost && parent2isGhost {
			gsToSelfHistMap[igs] = gsToSelfHistMap[histElem.parent1]
			continue
		}

		if histElem.parent2 >= 0 {
			gsToSelfHistMap[igs] = len(aa.cs.history)

			jeti := aa.cs.history[gsToSelfHistMap[histElem.parent1]].jet
			jetj := aa.cs.history[gsToSelfHistMap[histElem.parent2]].jet

			_, err := aa.cs.ijRecombinationStep(jeti, jetj, histElem.dij)
			if err != nil {
				return nil, err
			}

		} else {
			// check if parent 2 is beam

			gsToSelfHistMap[igs] = len(aa.cs.history)
			err := aa.cs.ibRecombinationStep(aa.cs.history[gsToSelfHistMap[histElem.parent1]].jet, histElem.dij)
			if err != nil {
				return nil, err
			}
		}
	}
	return gs, nil
}

func (aa *ClusterSequenceActiveArea) postProcess() {
	for _, a := range aa.averageArea {
		a /= float64(aa.ghostSpec.repeat)
	}
	for _, a2 := range aa.averageArea2 {
		a2 /= float64(aa.ghostSpec.repeat)
	}
	if aa.ghostSpec.repeat > 1 {
		temp := float64(aa.ghostSpec.repeat - 1)
		for i := range aa.averageArea2 {
			aa.averageArea2[i] = math.Sqrt(math.Abs(aa.averageArea2[i]-aa.averageArea[i]*aa.averageArea[i]) / temp)
		}

	} else {
		for i := range aa.averageArea2 {
			aa.averageArea2[i] = 0
		}
	}
	aa.nonJetArea /= float64(aa.ghostSpec.repeat)
	aa.nonJetArea2 /= float64(aa.ghostSpec.repeat)
	aa.nonJetArea2 = math.Sqrt(math.Abs(aa.nonJetArea2-aa.nonJetArea*aa.nonJetArea) / float64(aa.ghostSpec.repeat))
	aa.nonJetN /= float64(aa.ghostSpec.repeat)

	for _, a := range aa.averageP4 {
		x := 1 / float64(aa.ghostSpec.repeat) * a.Px()
		y := 1 / float64(aa.ghostSpec.repeat) * a.Py()
		z := 1 / float64(aa.ghostSpec.repeat) * a.Pz()
		e := 1 / float64(aa.ghostSpec.repeat) * a.E()
		p4 := fmom.NewPxPyPzE(x, y, z, e)
		a.Set(p4.Clone())
	}
}

type ClusterSequenceActiveAreaExplicitGhosts struct {
	cs           *ClusterSequence
	ghosts       []Jet // FIXME probably not needed. If needed something needs to be added
	ghostArea    float64
	jets         []Jet
	isPureGhost  []bool
	initialHardN int
	nGhosts      int
	areas        []float64
	area4vectors []Jet
	def          JetDefinition
}

func NewClusterSequenceActiveAreaExplicitGhosts(jets []Jet, def JetDefinition,
	ghostspec *GhostedAreaSpec) (*ClusterSequenceActiveAreaExplicitGhosts, error) {
	eg := ClusterSequenceActiveAreaExplicitGhosts{
		def: def,
	}
	eg.ghosts = make([]Jet, 0)
	eg.isPureGhost = make([]bool, 0)
	eg.jets = make([]Jet, 0)
	err := eg.init(jets, def, ghostspec, eg.ghosts, 0)
	if err != nil {
		return nil, err
	}
	return &eg, nil
}

func (eg *ClusterSequenceActiveAreaExplicitGhosts) init(jets []Jet, def JetDefinition, ghostspec *GhostedAreaSpec,
	ghosts []Jet, ghostArea float64) error {
	for _, jet := range jets {
		eg.jets = append(eg.jets, jet)
		eg.isPureGhost = append(eg.isPureGhost, false)
	}
	eg.initialHardN = len(jets)
	if ghostspec != nil {
		eg.addGhostsSpec(ghostspec)
	} else {
		eg.addGhosts(ghosts, ghostArea)
	}
	var err error
	eg.cs, err = NewClusterSequence(eg.jets, def)
	if err != nil {
		return err
	}
	eg.postProcess()
	return nil
}

func (eg *ClusterSequenceActiveAreaExplicitGhosts) addGhosts(ghosts []Jet, ghostArea float64) {
	for _, ghost := range ghosts {
		eg.isPureGhost = append(eg.isPureGhost, true)
		eg.jets = append(eg.jets, ghost)
	}
	eg.ghostArea = ghostArea
	eg.nGhosts = len(ghosts)
}

func (eg *ClusterSequenceActiveAreaExplicitGhosts) addGhostsSpec(ghostspec *GhostedAreaSpec) {
	eg.jets = ghostspec.addGhosts(eg.jets)
	for i := eg.initialHardN; i < len(eg.jets); i++ {
		eg.isPureGhost = append(eg.isPureGhost, true)
	}
	eg.ghostArea = ghostspec.actualGhostArea
	eg.nGhosts = ghostspec.nGhosts
}

func (eg *ClusterSequenceActiveAreaExplicitGhosts) postProcess() error {
	eg.areas = make([]float64, len(eg.cs.history))
	eg.area4vectors = make([]Jet, len(eg.cs.history))
	for i := 0; i < eg.cs.initn; i++ {
		if eg.isPureGhost[i] {
			eg.areas[i] = eg.ghostArea
			eg.area4vectors[i].Set(eg.jets[i].Clone())
			x := eg.area4vectors[i].Px() * eg.ghostArea / eg.jets[i].Px()
			y := eg.area4vectors[i].Py() * eg.ghostArea / eg.jets[i].Py()
			z := eg.area4vectors[i].Pz() * eg.ghostArea / eg.jets[i].Pz()
			e := eg.area4vectors[i].E() * eg.ghostArea / eg.jets[i].E()
			p4 := fmom.NewPxPyPzE(x, y, z, e)
			eg.area4vectors[i].Set(p4.Clone())
		} else {
			eg.areas[i] = 0
			eg.area4vectors[i] = NewJet(0, 0, 0, 0)
		}
	}

	for i := eg.cs.initn; i < len(eg.cs.history); i++ {
		if eg.cs.history[i].parent2 == beamJetIndex {
			eg.isPureGhost = append(eg.isPureGhost,
				eg.isPureGhost[eg.cs.history[i].parent1])
			eg.areas[i] = eg.areas[eg.cs.history[i].parent1]
			eg.area4vectors[i] = eg.area4vectors[eg.cs.history[i].parent1]
		} else {
			eg.isPureGhost = append(eg.isPureGhost,
				eg.isPureGhost[eg.cs.history[i].parent1] &&
					eg.isPureGhost[eg.cs.history[i].parent2])
			eg.areas[i] = eg.areas[eg.cs.history[i].parent1] + eg.areas[eg.cs.history[i].parent2]
			var err error
			eg.area4vectors[i], err = eg.def.recombiner.Recombine(&eg.area4vectors[eg.cs.history[i].parent1],
				&eg.area4vectors[eg.cs.history[i].parent2])
			if err != nil {
				return err
			}

		}
	}
	return nil
}
