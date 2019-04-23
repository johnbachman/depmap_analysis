import csv
import pickle
from indra.tools import assemble_corpus as ac
from indra_db.client import get_statements_by_gene_role_type
from indra_db.util import get_primary_db, get_raw_stmts_frm_db_list


def get_increase_amt_stmts(filename):
    inc_stmts = get_statements_by_gene_role_type(
                        stmt_type='IncreaseAmount', fix_refs=False,
                        preassembled=False,
                        with_evidence=True, with_support=False)
    with open(filename, 'wb') as f:
        pickle.dump(inc_stmts, f)
    return inc_stmts


def get_decrease_amt_stmts(filename):
    dec_stmts = get_statements_by_gene_role_type(
                        stmt_type='DecreaseAmount', fix_refs=False,
                        preassembled=False,
                        with_evidence=True, with_support=False)
    with open(filename, 'wb') as f:
        pickle.dump(dec_stmts, f)
    return dec_stmts

def filter_stmts_for_keywords(stmts):
    kwds = {
        'kd': ['siRNA', 'silencing', 'si-', 'Si-', 'sh-', 'SH-', 'shRNA'],
        'ko': ['knockout', '-/-', 'KO'],
        'dn': ['dn-', 'DN-', 'dominant negative', 'dominant-negative'],
        'oe': ['overexpress', 'OE'],
        'chem': ['chemical inhibition of', 'inhibitor of']}
    results = {}
    for expt_type, kwd_list in kwds.items():
        for kwd in kwd_list:
            kwd_stmts = []
            for stmt in stmts:
                text = stmt.evidence[0].text
                if text and kwd in text:
                    kwd_stmts.append(stmt)
            results[kwd] = kwd_stmts
    return results


def dump_csv(filt_stmts, csv_filename):
    header =['expt_keyword', 'indra_stmt', 'subj_hgnc', 'subj_text',
             'obj_hgnc', 'obj_text', 'sentence']
    rows = [header]
    for kwd, stmt_list in filt_stmts.items():
        for s in stmt_list:
            text = s.evidence[0].text
            if not (text and s.subj.db_refs.get('HGNC') and
                             s.obj.db_refs.get('HGNC')):
                continue
            row = [kwd, str(s),
                   s.subj.db_refs['HGNC'], s.subj.db_refs['TEXT'],
                   s.obj.db_refs['HGNC'], s.obj.db_refs['TEXT'],
                   text]
            rows.append(row)
    with open(csv_filename, 'wt') as f:
        csvwriter = csv.writer(f, delimiter=',')
        csvwriter.writerows(rows)


if __name__ == '__main__':
    """
    reload = False
    if reload:
        print("Getting IncreaseAmount statements")
        inc_stmts = get_increase_amt_stmts('increase_amt.pkl')
        print("Getting DecreaseAmount statements")
        dec_stmts = get_decrease_amt_stmts('decrease_amt.pkl')
        stmts = inc_stmts + dec_stmts
    else:
        inc_stmts = ac.load_statements('increase_amt.pkl')
        dec_stmts = ac.load_statements('decrease_amt.pkl')
        stmts = inc_stmts + dec_stmts
    """
    stmts = ac.load_statements('stmt_sample.pkl')
    stmts = [s for s in stmts if s.evidence[0].source_api == 'reach']
    stmts = [s for s in stmts if None not in s.agent_list()]
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_genes_only(stmts)
    stmts = ac.filter_human_only(stmts)
    filt_stmts = filter_stmts_for_keywords(stmts)
    dump_csv(filt_stmts, 'kwd_stmts.csv')
