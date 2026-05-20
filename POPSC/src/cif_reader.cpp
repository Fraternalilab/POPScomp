/*==============================================================================
cif_reader.cpp : C++ code to process mmcif structures 

Copyright (C) 2026 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/model.hpp>

#include <zlib.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <string>
#include <exception>

#include "cif_reader.h"

/*___________________________________________________________________________*/
static char* copy_string(const std::string& x) {
    char* p = (char*) std::malloc(x.size() + 1);
    if (!p) return nullptr;
    std::memcpy(p, x.c_str(), x.size() + 1);
    return p;
}

/*___________________________________________________________________________*/
static std::string gunzip_to_temp_file(const char* filename) {
    char tmpname[] = "/tmp/pops_mmcif_XXXXXX";

    int fd = mkstemp(tmpname);
    if (fd < 0)
        return "";

    gzFile in = gzopen(filename, "rb");
    if (!in) {
        close(fd);
        unlink(tmpname);
        return "";
    }

    char buffer[8192];
    int nread;

    while ((nread = gzread(in, buffer, sizeof(buffer))) > 0) {
        if (write(fd, buffer, nread) != nread) {
            gzclose(in);
            close(fd);
            unlink(tmpname);
            return "";
        }
    }

    gzclose(in);
    close(fd);

    return std::string(tmpname);
}

/*___________________________________________________________________________*/
Structure* read_cif(const char* filename) {
    Structure* s = nullptr;

    try {
        gemmi::cif::Document doc;
        std::string fname(filename);

        if (fname.size() >= 3 &&
            fname.substr(fname.size() - 3) == ".gz") {

            std::string tmpfile = gunzip_to_temp_file(filename);
            if (tmpfile.empty())
                return nullptr;

            doc = gemmi::cif::read_file(tmpfile);
            unlink(tmpfile.c_str());

        } else {
            doc = gemmi::cif::read_file(fname);
        }

        if (doc.blocks.empty())
            return nullptr;

        gemmi::Structure st =
            gemmi::make_structure_from_block(doc.blocks[0]);

        if (st.models.empty())
            return nullptr;

        const gemmi::Model& model = st.models[0];

        int n = 0;
        for (const gemmi::Chain& chain : model.chains) {
            for (const gemmi::Residue& res : chain.residues) {
                n += (int) res.atoms.size();
            }
        }

        s = (Structure*) std::calloc(1, sizeof(Structure));
        if (!s)
            return nullptr;

        s->natom = n;

        s->xyz         = (double*) std::malloc(sizeof(double) * 3 * n);
        s->atom_name   = (char**) std::calloc(n, sizeof(char*));
        s->atom_number = (int*)    std::malloc(sizeof(int) * n);
        s->altloc      = (char*)   std::malloc(sizeof(char) * n);
        s->res_name    = (char**) std::calloc(n, sizeof(char*));
        s->res_number  = (int*)    std::malloc(sizeof(int) * n);
        s->ins_code    = (char*)   std::malloc(sizeof(char) * n);
        s->chain_name  = (char**) std::calloc(n, sizeof(char*));
        s->element     = (char**) std::calloc(n, sizeof(char*));
        s->record_type = (char*)   std::malloc(sizeof(char) * n);

        if (!s->xyz || !s->atom_name || !s->atom_number ||
            !s->altloc || !s->res_name || !s->res_number ||
            !s->ins_code || !s->chain_name || !s->element ||
            !s->record_type) {
            free_structure(s);
            return nullptr;
        }

        int i = 0;
        s->chain_number = 0;
        s->nresidue = 0;

        for (const gemmi::Chain& chain : model.chains) {
            s->nresidue += (int) chain.residues.size();

            for (const gemmi::Residue& res : chain.residues) {
                for (const gemmi::Atom& atom : res.atoms) {

                    s->xyz[3*i + 0] = atom.pos.x;
                    s->xyz[3*i + 1] = atom.pos.y;
                    s->xyz[3*i + 2] = atom.pos.z;

                    s->atom_name[i] = copy_string(atom.name);
                    s->atom_number[i] = atom.serial;
                    s->altloc[i] = atom.altloc ? atom.altloc : ' ';

                    s->res_name[i] = copy_string(res.name);
                    s->res_number[i] =
                        res.seqid.num.has_value() ? res.seqid.num.value : 0;
                    s->ins_code[i] = res.seqid.icode ? res.seqid.icode : ' ';

                    s->chain_name[i] = copy_string(chain.name);

                    std::string elem = atom.element.name();
                    if (elem.empty() || elem == "X")
                        elem = atom.name;

                    s->element[i] = copy_string(elem);

                    s->record_type[i] = res.het_flag;

                    if (!s->atom_name[i] || !s->res_name[i] ||
                        !s->chain_name[i] || !s->element[i]) {
                        free_structure(s);
                        return nullptr;
                    }

                    ++i;
                }
            }

            ++s->chain_number;
        }

        return s;
    }

    catch (const std::exception& e) {
        fprintf(stderr, "Gemmi CIF read error: %s\n", e.what());
        if (s)
            free_structure(s);
        return nullptr;
    }
}

/*___________________________________________________________________________*/
void free_structure(Structure* s) {
    if (!s)
        return;

    if (s->atom_name) {
        for (int i = 0; i < s->natom; ++i)
            std::free(s->atom_name[i]);
    }

    if (s->res_name) {
        for (int i = 0; i < s->natom; ++i)
            std::free(s->res_name[i]);
    }

    if (s->chain_name) {
        for (int i = 0; i < s->natom; ++i)
            std::free(s->chain_name[i]);
    }

    if (s->element) {
        for (int i = 0; i < s->natom; ++i)
            std::free(s->element[i]);
    }

    std::free(s->xyz);
    std::free(s->atom_name);
    std::free(s->atom_number);
    std::free(s->altloc);
    std::free(s->res_name);
    std::free(s->res_number);
    std::free(s->ins_code);
    std::free(s->chain_name);
    std::free(s->element);
    std::free(s->record_type);
    std::free(s);
}

/*============================================================================*/

