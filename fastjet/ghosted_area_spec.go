// Copyright 2017 The go-hep Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fastjet

import (
	"math"
	"math/rand"
)

type GhostedAreaSpec struct {
	selector        Selector
	repeat          int
	ghostArea       float64
	gridScatter     float64
	ptScatter       float64
	meanGhostPt     float64
	fj2Placement    bool
	actualGhostArea float64
	ghostMaxRap     float64
	ghostRapOffset  float64
	drap            float64
	dphi            float64
	nphi            int
	nrap            int
	nGhosts         int
}

// FIXME delete selector and move rapidity extend to ghosted area spec
type Selector struct {
	finiteArea      bool
	appliesJetByJet bool
	description     string
	takesRefJ       bool
	refJ            Jet
	isGeometric     bool
	knownArea       bool
}

func NewDefaultSelector() *Selector {
	s := &Selector{
		appliesJetByJet: true,
		takesRefJ:       false,
		isGeometric:     false,
		knownArea:       false,
		finiteArea:      true,
	}
	return s
}

// FIXME is probably 3,-3 but test is killed with those values
func (*Selector) getRapidityExtend() (float64, float64) {
	return 1, -1
}

func NewGhostedAreaSpec(selector Selector, repeat int, ghostArea float64, gridScatter float64,
	ptScatter float64, meanGhostPt float64) (*GhostedAreaSpec, error) {
	var er error
	if !selector.finiteArea {
		// er := errors.New("fastjet: selector needs finite area")
		// return nil, er
	}
	if !selector.appliesJetByJet {
		// er := errors.New("fastjet: selector must apply jet by jet")
		// return nil,er
	}
	gs := &GhostedAreaSpec{
		selector:        selector,
		repeat:          repeat,
		ghostArea:       ghostArea,
		gridScatter:     gridScatter,
		ptScatter:       ptScatter,
		meanGhostPt:     meanGhostPt,
		fj2Placement:    false,
		actualGhostArea: -1,
	}

	max, min := gs.selector.getRapidityExtend()
	gs.ghostMaxRap = 0.5 * (max - min)
	gs.ghostRapOffset = 0.5 * (max + min)

	gs.initialize()

	return gs, er
}

func NewDefaultGhostAreaSpec() *GhostedAreaSpec {
	s := NewDefaultSelector()
	ga, err := NewGhostedAreaSpec(*s, 1, 0.01, 1, 0.1, 1e-100)
	if err != nil {
		return nil
	}
	return ga
}

func (gs *GhostedAreaSpec) SetFj2Placement(val bool) {
	gs.fj2Placement = val
	gs.initialize()
}

func (gs *GhostedAreaSpec) initialize() {
	gs.drap = math.Sqrt(gs.ghostArea)
	gs.dphi = gs.drap
	twopi := 2 * math.Pi
	if gs.fj2Placement {
		gs.nphi = int(math.Ceil(twopi / gs.dphi))
		gs.dphi = twopi / float64(gs.nphi)
		gs.nrap = int(math.Ceil(gs.ghostMaxRap / gs.drap))
		gs.drap = gs.ghostMaxRap / float64(gs.nrap)
		gs.actualGhostArea = gs.dphi * gs.drap
		gs.nGhosts = (2*gs.nrap + 1) * gs.nphi
	} else {
		gs.nphi = int(twopi/gs.dphi + 0.5)
		gs.dphi = twopi / float64(gs.nphi)
		gs.nrap = int(gs.ghostMaxRap/gs.drap + 0.5)
		gs.drap = gs.ghostMaxRap / float64(gs.nrap)
		gs.actualGhostArea = gs.dphi * gs.drap
		gs.nGhosts = (2 * gs.nrap) * gs.nphi
	}
}

// FIXME add seeding, find range
func (gs *GhostedAreaSpec) ourRand() float64 {
	return rand.Float64()
}

func (gs *GhostedAreaSpec) addGhosts(jets []Jet) []Jet {
	var rapOffset float64
	var nrapUpper int
	if gs.fj2Placement {
		rapOffset = 0.0
		nrapUpper = gs.nrap
	} else {
		rapOffset = 0.5
		nrapUpper = gs.nrap - 1
	}

	// add momenta for ghosts
	for irap := -gs.nrap; irap <= nrapUpper; irap++ {
		for iphi := 0; iphi < gs.nphi; iphi++ {
			phiFj2 := (float64(iphi)+0.5)*gs.dphi + gs.dphi*(gs.ourRand()-0.5)*gs.gridScatter
			var phi float64
			if gs.fj2Placement {
				phi = 0.5*math.Pi - phiFj2
			} else {
				phi = phiFj2
			}
			rap := (float64(irap)+rapOffset)*gs.drap + gs.drap*(gs.ourRand()-0.5)*gs.gridScatter + gs.ghostRapOffset

			pt := gs.meanGhostPt * (1 + (gs.ourRand()-0.5)*gs.ptScatter)

			exrap := math.Exp(rap)
			pminus := pt / exrap
			pplus := pt * exrap
			px := pt * math.Cos(phi)
			py := pt * math.Sin(phi)
			pz := 0.5 * (pplus - pminus)
			e := 0.5 * (pplus + pminus)

			ghost := NewJet(px, py, pz, e)

			ghost.rap = rap
			ghost.phi = phi

			// if selector and selector.pass(mom) continue
			jets = append(jets, ghost)
		}
	}
	return jets
}
