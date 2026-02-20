/*==============================================================================
cif_reader.cpp : C++ code to process mmcif structures 
Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <gemmi/cif.hpp>     // gemmi::cif::Document, gemmi::cif::read_file
#include <gemmi/mmcif.hpp>   // gemmi::make_structure_from_block
#include <gemmi/model.hpp>   // gemmi::Structure, gemmi::Model, cra.atom
#include <cstdlib>
#include <string>

#include "cif_reader.h"

Structure* read_cif(const char* filename) {
    // IMPORTANT: use read_file(), NOT read(filename)
    gemmi::cif::Document doc = gemmi::cif::read_file(std::string(filename));

    if (doc.blocks.empty())
        return nullptr;

    // Convert mmCIF block -> Structure
    gemmi::Structure st = gemmi::make_structure_from_block(doc.blocks[0]);

    if (st.models.empty())
        return nullptr;

    const gemmi::Model& model = st.models[0];

    // Count atoms portably (older gemmi may not have Model::count_atoms())
    int n = 0;
	for (auto it = model.all().begin(); it != model.all().end(); ++it)
    	++n;

    Structure* s = (Structure*) std::malloc(sizeof(Structure));
    if (!s) return nullptr;

    s->natom = n;
    s->xyz = (double*) std::malloc(sizeof(double) * 3 * n);
    if (!s->xyz) {
        std::free(s);
        return nullptr;
    }

    int i = 0;
    for (auto cra : model.all()) {
        // In your gemmi, cra.atom is a pointer (const gemmi::Atom*)
        const gemmi::Atom* atom = cra.atom;
        s->xyz[3*i + 0] = atom->pos.x;
        s->xyz[3*i + 1] = atom->pos.y;
        s->xyz[3*i + 2] = atom->pos.z;
        ++i;
    }

    return s;
}

void free_structure(Structure* s) {
    if (!s) return;
    std::free(s->xyz);
    std::free(s);
}


/*============================================================================*/


