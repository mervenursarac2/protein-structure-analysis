import py3Dmol

with open("../data/1REX.pdb") as f:
    pdb_text = f.read()

view = py3Dmol.view(width=800, height=600)
view.addModel(pdb_text, "pdb")

# Genel yapıyı ikincil yapı tipine göre renklendir (heliks/tabaka/loop ayrı renk)
view.setStyle({"cartoon": {"colorscheme": "chainHetatm"}})
view.setStyle({}, {"cartoon": {"color": "lightgray"}})  # önce hepsini gri yap

# Heliksleri ve tabakaları otomatik tanıyıp renklendir
view.setStyle({}, {"cartoon": {"color": "white"}})

# Pozisyon 56'yı (mutasyon bölgesi) KIRMIZI ve "stick" (atom detaylı) göster
view.addStyle({"resi": "56"}, {"stick": {"color": "red"}, "sphere": {"scale": 0.3, "color": "red"}})

# Komşu bölgeyi (54-58) sarı ile vurgula, bağlam için
view.addStyle({"resi": "54-58"}, {"cartoon": {"color": "yellow"}})

view.zoomTo({"resi": "56"})  # kamerayı mutasyon bölgesine odakla
view.write_html("../results/lysozyme_I56_highlight.html")
print("Kaydedildi: results/lysozyme_I56_highlight.html")