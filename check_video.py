import json

with open('notebooks/party_cartograms_zweitstimme_prozent.ipynb') as f:
    nb = json.load(f)
    for i, cell in enumerate(nb['cells']):
        src = ''.join(cell.get('source', []))
        outputs = cell.get('outputs', [])

        # Check source
        if 'video' in src.lower() or 'mp4' in src.lower():
            print(f'Cell {i} SOURCE:')
            print(src[:500])
            print('---')

        # Check outputs
        for out in outputs:
            out_str = str(out)
            if 'video' in out_str.lower() or 'mp4' in out_str.lower():
                print(f'Cell {i} OUTPUT:')
                print(out_str[:500])
                print('---')
