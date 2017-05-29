// Copyright 2017 The go-hep Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fastjet

import (
	"errors"
)

// FIXME original Problem: all areas are 0.
// it looks like the ghosts are not colliding with the jets and since
// only the ghost have areas and then the parents of the ghosts, none
// of the jets have areas

// after changing up some values in the ghostedareaspec some jets have areas now,
// but the values are still incorrect

type ClusterSequenceArea struct {
	cs   *ClusterSequence
	area AreaDefinition
	aa   *ClusterSequenceActiveArea
}

func NewClusterSequenceArea(jets []Jet, def JetDefinition, area AreaDefinition) (*ClusterSequenceArea, error) {
	//not needed here
	cs, err := NewClusterSequence(jets, def)
	if err != nil {
		return nil, err
	}

	csa := &ClusterSequenceArea{
		cs:   cs,
		area: area,
	}
	if area.AreaType == ActiveArea {
		aa, err := NewClusterSequenceActiveArea(jets, def, area.GhostSpec)
		if err != nil {
			return nil, err
		}
		csa.aa = aa
		err = csa.aa.initializeAndRunAA()
		if err != nil {
			return nil, err
		}
	} else {
		err := errors.New("fastjet: area type not implemented")
		return nil, err
	}
	return csa, nil
}

// FIXME make it all general so it works for all area types
func (csa *ClusterSequenceArea) Area(jet *Jet) float64 {
	return csa.aa.averageArea[jet.hidx]
}

func (csa *ClusterSequenceArea) AreaErr(jet *Jet) float64 {
	return csa.aa.averageArea2[jet.hidx]
}

func (csa *ClusterSequenceArea) NumExclusiveJets(dcut float64) int {
	return csa.cs.NumExclusiveJets(dcut)
}

func (csa *ClusterSequenceArea) ExclusiveJets(dcut float64) ([]Jet, error) {
	return csa.cs.ExclusiveJets(dcut)
}

func (csa *ClusterSequenceArea) ExclusiveJetsUpTo(njets int) ([]Jet, error) {
	return csa.cs.ExclusiveJetsUpTo(njets)
}

func (csa *ClusterSequenceArea) InclusiveJets(ptmin float64) ([]Jet, error) {
	return csa.cs.InclusiveJets(ptmin)
}
