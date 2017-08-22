// Copyright 2017 The go-hep Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package heap

import (
	"testing"
)

func TestHeapPush(t *testing.T) {
	h := NewHeap()
	h.Push(0, 1, 10)
	h.Push(2, 3, 5)
	h.Push(4, 5, 8)
	h.Push(6, 7, 2)
	want := []struct {
		jeti, jetj int
	}{
		{-1, -1},
		{6, 7},
		{2, 3},
		{4, 5},
		{0, 1},
	}
	got := h.items
	if len(got) != len(want) {
		t.Errorf("got = %d items, want = %d", len(got), len(want))
	}
	for i, item := range got {
		if i == 0 {
			continue
		}
		if item.jeti != want[i].jeti {
			t.Errorf("h.items[%d].jeti: got = %d, want = %d", i, item.jeti, want[i].jeti)
		}
		if item.jetj != want[i].jetj {
			t.Errorf("h.items[%d].jetj: got = %d, want = %d", i, item.jetj, want[i].jetj)
		}
	}
}

func TestHeapPop(t *testing.T) {
	items := []item{item{}, item{jeti: 6, jetj: 7, dij: 2}, item{jeti: 2, jetj: 3, dij: 5}, item{jeti: 4, jetj: 5, dij: 8},
		item{jeti: 0, jetj: 1, dij: 10}, item{jeti: 8, jetj: 9, dij: 6},
	}
	h := &Heap{items: items, n: 5}
	want := []struct {
		jeti, jetj int
	}{
		{6, 7},
		{2, 3},
		{8, 9},
		{4, 5},
		{0, 1},
	}
	got := make([]struct{ jeti, jetj int }, len(items)-1)
	for i := 0; !h.IsEmpty(); i++ {
		if i >= len(got) {
			t.Fatalf("Heap with n = %d should be empty", h.n)
		}
		jeti, jetj := h.Pop()
		got[i].jeti = jeti
		got[i].jetj = jetj
	}
	if len(got) != len(want) {
		t.Errorf("got = %d items, want = %d", len(got), len(want))
	}
	for i := range got {
		if want[i].jeti != got[i].jeti {
			t.Errorf("got[%d].jeti = %d, want[%d].jeti = %d", i, got[i].jeti, i, want[i].jeti)
		}
		if want[i].jetj != got[i].jetj {
			t.Errorf("got[%d].jetj = %d, want[%d].jetj = %d", i, got[i].jetj, i, want[i].jetj)
		}
	}
}
