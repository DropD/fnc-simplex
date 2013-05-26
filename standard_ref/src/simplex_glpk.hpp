#include <glpk.h>

template <typename T>
class SimplexGLPK {
    private:
    glp_prob *lp;
    glp_smcp parm;

    public:
    unsigned int PERFC_MEM, PERFC_ADDMUL, PERFC_DIV;

    SimplexGLPK() {
        glp_term_out(GLP_OFF);
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_OFF;
    }
    void load(std::string fname) {
        lp = glp_create_prob();
        glp_read_lp(lp, NULL, fname.c_str());
    }

    inline void solve() {
        glp_simplex(lp, &parm);
    }
};
