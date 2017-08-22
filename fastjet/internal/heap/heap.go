// Copyright 2017 The go-hep Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package heap

// item is an item in the heap. It contains the indices of the two jets to be clustered
// and their distance in the Eta-Phi plane.
type item struct {
	// jeti and jetj are the indices of the jets to be recombined.
	jeti, jetj int
	// dij is the distance between the jets in the Eta-Phi plane.
	dij float64
}

// Heap contains a slice of items and the last index.
type Heap struct {
	n     int // n is the last index of the slice
	items []item
}

// NewHeap returns a heap pointer.
func NewHeap() *Heap {
	h := &Heap{n: 0}
	// item{} is a placeholder at h.items[0]
	h.items = append(h.items, item{})
	return h
}

// Push inserts two new clustering candidates and their distance.
func (h *Heap) Push(jeti, jetj int, dij float64) {
	item := item{jeti: jeti, jetj: jetj, dij: dij}
	h.n++
	if h.n >= len(h.items) {
		h.items = append(h.items, item)
	} else {
		h.items[h.n] = item
	}
	h.moveUp(h.n)
}

// Pop returns the two jets with the lowest distance.
// It returns -1, -1 when the heap is empty.
func (h *Heap) Pop() (jeti, jetj int, dij float64) {
	if h.n == 0 {
		return -1, -1
	}
	item := h.items[1]
	h.n--
	if h.n == 0 {
		return item.jeti, item.jetj
	}
	h.swap(1, h.n+1)
	h.moveDown(1)
	return item.jeti, item.jetj, item.dij
}

// IsEmpty returns whether a heap is empty.
func (h *Heap) IsEmpty() bool {
	return h.n == 0
}

// moveUp compares an item with its parent and moves it up if it has
// a lower distance. It will keep moving up until the parent's distance is less
// or it reaches the top.
func (h *Heap) moveUp(i int) {
	parent := int(i / 2)
	if parent == 0 {
		return
	}
	if h.less(i, parent) {
		h.swap(i, parent)
		h.moveUp(parent)
	}
}

// moveDown compares an item to its children and moves it down
// if one of the children has a lower distance. It will keep
// moving down until it reaches a leaf or both children have
// a higher distance.
func (h *Heap) moveDown(i int) {
	leftChid := 2 * i
	if leftChid > h.n {
		return
	}
	rightChild := 2*i + 1
	var smallestChild int
	if leftChid == h.n {
		smallestChild = leftChid
	} else if h.less(leftChid, rightChild) {
		smallestChild = leftChid
	} else {
		smallestChild = rightChild
	}
	if h.less(smallestChild, i) {
		h.swap(i, smallestChild)
		h.moveDown(smallestChild)
	}
}

// less returns true if the item at index i has a lower distance
// than the item at index j.
func (h *Heap) less(i, j int) bool {
	return h.items[i].dij < h.items[j].dij
}

// swap swaps the two items at the indices i and j.
func (h *Heap) swap(i, j int) {
	h.items[i], h.items[j] = h.items[j], h.items[i]
}
