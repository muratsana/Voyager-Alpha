# Voyager Alpha — GUI Yeniden Tasarım Spesifikasyonu

**Amaç:** Mevcut PyQt6 arayüzünü (`gui/asteroid_workspace.py` 1557 satır, `gui/exoplanet_workspace.py` 1186 satır, `gui/theme.py`) kullanıcı için anlaşılır, bilimsel iş akışını yönlendiren ve raporlamaya kadar götüren bir yapıya dönüştürmek. Bu doküman **tasarım kararlarını, ekran yapısını, bileşen sözleşmelerini ve QSS/token tablosunu** içerir; AI coder doğrudan uygulayabilir.

---

## 1. Mevcut arayüzün sorunları (kod incelemesinden)

| # | Sorun | Kanıt | Etki |
|---|---|---|---|
| G1 | **Her şey tek ekranda, 13+ bölüm (section) aynı anda görünür.** Sol panelde Sequence Validation + Workflow Rail + Calibration + Detection Profile + Manual Commands; sağda Evidence + Cross-match + Review + Confirmation & Export; ortada viewer + filmstrip + blink + log + sonuç tabloları. | `asteroid_workspace._build_left_panel/_build_right_panel/_build_center_panel` | Kullanıcı nereden başlayacağını bilmiyor; adımlar arası bağımlılık (önce klasör, sonra doğrula, sonra WCS…) görünmüyor. |
| G2 | **Sabit piksel genişlikleri ve yükseklikler.** Sol panel 260–320 px, sağ 300–390 px, sonuç bölmesi `setFixedHeight(230)`, notlar 42 px, butonlar 102 px, min pencere 1280×760. | `setMinimumWidth/MaximumWidth/setFixedHeight` çağrıları | 1366×768 dizüstülerde görüntü alanı ~700×300 px'e düşüyor; 4K'da paneller minicik kalıyor; DPI ölçekleme yok. |
| G3 | **Tipografi çok küçük ve tutarsız.** Global 12 px, caption 9 px, tablo başlığı 10 px, badge 10 px; "Segoe UI" sabit (Windows dışı yok). | `theme.py` | Okunabilirlik düşük; WCAG kontrast bazı renklerde (#78909a / #10171b ≈ 3.6:1) yetersiz. |
| G4 | **Dil karışık:** "Sequence Validation", "Kaynak Klasör", "Guided Analizi Başlat", "Discover Unknown Movers", "Seçili Objeyi Ortala", "Needs follow-up". | tüm GUI | Profesyonel görünmüyor; çeviri altyapısı yok. |
| G5 | **Modal QMessageBox ile ayar/yardım.** Ayarlar penceresi yalnızca "ASTAP bulundu/bulunamadı" mesajı; gözlemevi/kamera/teleskop profili yok. | `main_window.show_settings` | Zorunlu bilimsel girdiler (konum, MPC kodu, doğrusallık) girilecek yer yok. |
| G6 | **Parametreler kullanıcıya anlamsız birimlerle sunuluyor:** "1.5 px seed hareketi", "2.8 px yeniden eşleme", "0.9/1.8 px RMS", "Aperture radius 6 px", "sigma 5". | `DetectionSettingsDialog`, `combo_aperture` | Kullanıcı ″/dk, FWHM katı gibi fiziksel birim bekler; profil (MBA/NEO/TNO) yok. |
| G7 | **Durum ve hata iletişimi log tablosuna gömülü.** `INFO|WARN|ERROR` satırları; hangi adımın neden başarısız olduğu, ne yapılması gerektiği yok. | `AnalysisLog` | Kullanıcı log okumak zorunda. |
| G8 | **Işık eğrisi tek panel** (Relative flux vs Elapsed time); artık/binli/karşılaştırma/sistematik yok; asteroid tarafında hız–PA grafiği, RMS-vs-bin yok. | `exoplanet_workspace.plot_result` | Bilimsel değerlendirme yapılamıyor. |
| G9 | **Overlay etkileşimi zayıf:** tıklanan yıldızın bilgisi yok, katman yöneticisi 4 checkbox, hareket vektörü yok, RA/Dec–piksel okuma yok. | `FitsViewer` | Kullanıcı isteği (bilinen nesneleri overlay'de, hareket eksenini göster) karşılanmıyor. |
| G10 | **Durum kalıcılığı yok:** son klasör, profil, pencere düzeni, seçilen yıldızlar kaydedilmiyor (`QSettings` yok). | — | Her açılışta sıfırdan. |
| G11 | **Klavye kısayolu, undo, sürükle-bırak klasör, çoklu dizi yok.** | — | Verimsiz. |

---

## 2. Tasarım ilkeleri

1. **İş akışı = arayüz.** Her modül, soldan sağa ilerleyen **adım şeridi** (stepper) ile yönetilir; yalnızca aktif adımın kontrolleri görünür, tamamlananlar özet satırına çöker (progressive disclosure).
2. **Görüntü merkezde, her zaman büyük.** Görüntü alanı pencerenin ≥ %55'i; yan paneller **dock** (kapatılabilir, yüzdürülebilir, kaydedilebilir).
3. **Fiziksel birimler.** Her sayısal kontrol birimli (″/dk, mag, s, FWHM×); dahili piksel değerleri türetilir ve yanında gri gösterilir.
4. **Kararı kullanıcı verir, program kanıt sunar.** Her sonuç kartı: "Ne bulundu / Ne kadar güvenilir / Neden / Sıradaki adım".
5. **Tek dil, çevrilebilir.** Tüm metinler `tr()` ile; varsayılan Türkçe, İngilizce çeviri dosyası (`.ts`). Bilimsel terimler parantez içinde İngilizce (ör. "Karşılaştırma yıldızı (comparison)").
6. **DPI ve pencere bağımsız.** Sabit piksel yok; `QSizePolicy`, oran, `em` tabanlı boşluk; min pencere 1200×700, hedef 1920×1080, 4K'da otomatik ölçek (`Qt.HighDpiScaleFactorRoundingPolicy.PassThrough`).
7. **Erişilebilirlik.** Metin/arka plan kontrastı ≥ 4.5:1; odak halkası; tüm eylemler klavyeyle; tooltip + durum çubuğu yardımı.

---

## 3. Uygulama iskeleti

```
MainWindow (QMainWindow)
├── Üst şerit (48 px): Logo · Modül sekmeleri [Asteroid | Ötegezegen] · Proje adı · Profil seçici (Gözlemevi/Kamera/Teleskop) · Ağ durumu (çevrimiçi/önbellek) · Ayarlar · Yardım
├── Adım şeridi (Stepper, 40 px): modüle özgü adımlar (aşağıda) — tıklanabilir, durum renkli (bekliyor/aktif/tamam/uyarı/hata)
├── Merkez: QMainWindow dock sistemi
│   ├── Merkez widget: ImageViewer (pan/zoom/crosshair/stretch/overlay katmanları/blink)
│   ├── Sol dock "Adım paneli": yalnızca aktif adımın form/kontrolleri (QStackedWidget)
│   ├── Sağ dock "Sonuç & Kanıt": seçili nesnenin kartı (tracklet / yıldız / ışık eğrisi noktası)
│   ├── Alt dock "Grafikler": sekmeli (Işık eğrisi | Artıklar | Karşılaştırmalar | Sistematikler | RMS-vs-bin) veya (Hız–PA | ξ/η fit | SNR-zaman | Yetenek kartı)
│   ├── Alt dock "Tablolar": Tracklet'ler / Bilinen nesneler / Kare ölçümleri / Kare kalitesi
│   └── Alt dock "Günlük": filtrelenebilir (INFO/WARN/ERROR), adım etiketli, "kopyala" düğmesi
└── Durum çubuğu (28 px): aşama ilerlemesi (segmentli) · kare i/N · imleç (x, y, ADU, RA/Dec) · zoom · BJD/UTC saati · ağ
```

Dock yerleşimi `QSettings` ile kaydedilir; "Yerleşimi sıfırla" menüsü.

### 3.1 Asteroid modülü adımları

| Adım | Panel içeriği | Tamamlanma koşulu | Özet satırı |
|---|---|---|---|
| 1. Gözlem | Klasör seç (sürükle-bırak), dosya listesi (zaman sıralı, kalite bayrakları), **Dizi yetenek kartı** (asteroid dok. §2.7), rejim seçici (MBA/NEO/TNO/Özel → hız penceresi ″/dk, PA aralığı) | ≥ 3 geçerli kare, zaman monoton | "24 kare · 62 dk · 0.78″/px · MBA rejimi" |
| 2. Kalibrasyon | Bias/Dark/Flat: ham klasörden **master üret** veya mevcut master seç; eşleşme kontrolleri (binning, ΔT, poz, filtre); kötü piksel maskesi önizleme | (isteğe bağlı) | "B+D+F · 17/20/25 kare · 143 kötü piksel" |
| 3. Astrometri | Her kare plate solve (ASTAP/astrometry.net), sonuç tablosu (RMS ″, yıldız sayısı), kabul eşiği kaydırıcı, başarısızlar için yeniden dene | ≥ 3 kare RMS ≤ 0.5″ | "24/24 çözüldü · medyan RMS 0.31″" |
| 4. Bilinen nesneler | SkyBoT/sb_ident sorgu (topo), liste (isim, V, hız, PA, Err), overlay'de beklenen yol okları, "görünür mü" açıklık SNR'ı | sorgu tamamlandı veya çevrimdışı atlandı | "17 bilinen · 11 görünür · 2 NEO" |
| 5. Tespit & Bağlama | Eşik (σ), sabit kaynak yöntemi (Gaia maskesi/medyan/fark), hız penceresi (rejimden), tutarlılık toleransları — **Basit/Gelişmiş** anahtarı | çalıştırıldı | "312 tespit · 9 tracklet · 2 bilinmeyen aday" |
| 6. İnceleme | Tracklet listesi + kanıt kartı + nesne-merkezli blink + sentetik izleme; Kabul/Ret/Takip; uydu (TLE) ve NEOCP kontrol düğmeleri; Find_Orb Väisälä RMS | her aday incelendi | "2 kabul · 5 ret · 2 takip" |
| 7. Rapor | ADES (PSV/XML, `submit.xsd` doğrulama), MPC80, CSV, HTML, ALCDEF; gönderim öncesi kontrol listesi (kaide dok. §12) | dosya yazıldı | "ADES: 6 gözlem · 2 tracklet" |

### 3.2 Ötegezegen modülü adımları

| Adım | Panel | Koşul |
|---|---|---|
| 1. Gözlem | Klasör, kare listesi, **hedef seçimi katalogdan** (NEA/ExoClock arama kutusu) veya tıklayarak; efemeris: T1/T4 ± σ, taban çizgisi yeterliliği, ExoClock önceliği; "çekim uygunluk kartı" (kadans, verim, tepe ADU, FWHM px, X aralığı) | ≥ 5 kare |
| 2. Kalibrasyon | (asteroid ile aynı bileşen) | — |
| 3. Astrometri | Plate solve; **tüm alan overlay**: bilinen/aday/FP ötegezegen ev sahipleri (katman renkleri), Gaia değişkenler, VSX | — |
| 4. Yıldız seçimi | Hedef T1; otomatik karşılaştırma önerisi (puan tablosu: Δmag, BP−RP, RMS, VSX, mesafe, ret nedeni); elle ekle/çıkar; T2… NEB yıldızları (SG1 modu) | ≥ 3 karşılaştırma (uyarı ile 1–2) |
| 5. Fotometri | Açıklık modu (otomatik FWHM taraması / sabit / değişken), halka, gökyüzü yöntemi; hata modeli (gain, RON, scintillation) ; çalıştır → kare ölçüm tablosu | — |
| 6. Model & Kalite | Detrend seçimi (AIRMASS önce, ΔBIC kuralı), fit (batman/pylightcurve + LD tablosu), MCMC/nested; kalite kartı: derinlik ± σ, T_mid ± σ (BJD_TDB), Rp/R*, RMS, β, χ²_red, ExoClock kabul testleri; grafik seti | — |
| 7. Rapor | AAVSO Exoplanet DB dosyası, ExoClock 3-sütun, ETD Dmag, PNG/PDF grafik, ölçüm tablosu (AIJ sütunları) | — |

---

## 4. Bileşen sözleşmeleri

### 4.1 ImageViewer
- Pan (orta tuş/boşluk+sürükle), zoom (tekerlek, imleç merkezli, 1:1 / sığdır / 200%), mini harita.
- İmleç okuma: x, y, ADU, RA/Dec (WCS varsa), en yakın katalog nesnesi.
- Stretch: Auto STF, asinh, log, lineer, manuel (histogram widget'ı ile), invert; **kareye özgü** ve **diziye kilitli** seçenekleri.
- Overlay katman yöneticisi (ağaç): Bilinen asteroid/kuyruklu yıldız (sınıf renkleri, belirsizlik elipsi, **hareket vektörü + kare noktaları**), tracklet'ler (durum renkleri, fit çizgisi, tolerans bandı), tespitler, ötegezegen ev sahipleri (durum renkleri), Gaia değişkenler, VSX, NEB 2.5′ dairesi, fotometri açıklıkları (T1/C2…), uydu geçişleri, ızgara, Kuzey/Doğu oku, ölçek çubuğu, maskeler. Her katman: görünürlük, opaklık, etiket aç/kapa.
- Tıklama modları: Seç (nesne kartı açılır) / Hedef seç / Karşılaştırma ekle / Ölç (geçici açıklık) / Mesafe ölç.
- Blink: kare listesi, fps, ileri/geri, **nesne-merkezli** (seçili tracklet veya bilinen nesne hızıyla kaydırarak), fark modu, döngü aralığı; klavye: ←/→ kare, boşluk oynat, B blink, D fark.
- Performans: piramit önbellek (QImage), thread'de stretch; 60 MB kare için ≤ 100 ms kare geçişi.

### 4.2 Sonuç kartı (sağ dock)
Şablon: Başlık (isim/ID) · Durum rozeti (Bilinen/Bilinmeyen aday/Ret/Uydu?) · **Ölçüler** (fiziksel birimlerle) · **Neden** (kural bazlı gerekçe listesi: "3 kare, RMS 0.4″, hız 0.62″/dk, SkyBoT eşleşmesi yok (en yakın 2024 AB, 48″)") · **Sıradaki adım** (öneri düğmeleri: "Sentetik izleme", "TLE kontrol", "NEOCP'de ara", "İkinci gece planla").

### 4.3 Parametre kontrolleri
- `UnitSpinBox`: değer + birim + türetilmiş piksel değeri gri; tooltip'te kaide referansı ("MOPS: ≤ 5°/gün").
- Basit/Gelişmiş anahtarı: Basit = rejim ön ayarı + 2–3 kontrol; Gelişmiş = tüm toleranslar.
- Her ayar değişikliğinde "Yeniden çalıştır" düğmesi kirli (dirty) durumuna geçer.

### 4.4 Grafik paneli (pyqtgraph)
Ötegezegen: (a) ışık eğrisi: ham gri nokta + 5 dk binli mavi + model kırmızı + T1/T4 kesikli + meridyen; (b) artıklar; (c) karşılaştırma eğrileri (offset); (d) sistematikler (FWHM, sky, airmass ters, x, y); (e) RMS-vs-bin log-log + 1/√N; (f) corner (MCMC). Tümü PNG/PDF dışa aktarılabilir, dok. §6.1 anotasyon kutusu ile.
Asteroid: hız–PA saçılımı (tespit edilenler + bilinen beklenenler), ξ(t)/η(t) fit + artık, SNR/kadir–zaman, plate-solve artık vektör alanı, shift-and-stack olabilirlik haritası, O−C tablosu.

### 4.5 Profil yöneticisi (Ayarlar)
Sekmeler: Gözlemevi (ad, enlem/boylam/yükseklik, MPC kodu, AAVSO kodu, saat dilimi, NTP/GPS notu), Teleskop (açıklık, odak, f/, plate scale doğrulama), Kamera (piksel, gain e⁻/ADU, RON, dark, full well, **doğrusallık sihirbazı**, Bayer), Filtreler (yerel ad → AAVSO/ADES band kodu), Araçlar (ASTAP yolu, astrometry.net indeks/API anahtarı, Find_Orb yolu), Ağ (önbellek yaşı, çevrimdışı mod), Görünüm (tema, dil, ölçek).

### 4.6 Boş durumlar ve hata metinleri
- Boş görüntü alanı: "Bir gözlem klasörü sürükleyin veya **Klasör Seç**. Desteklenen: FITS/FIT/FTS (2-B, 16/32-bit)."
- Hata kartı şablonu: **Ne oldu** · **Neden olabilir** · **Ne yapmalı** · "Günlüğü kopyala". Örn.: "Plate solve 3 karede başarısız — FOV tahmini yok (XPIXSZ/FOCALLEN eksik) → Teleskop profilinde plate scale girin veya `-fov` değerini ayarlayın."
- Uyarılar bloke etmez; RET seviyesindekiler adımı kırmızıya boyar ve ilerlemeyi kilitler (kaide dok. kontrol listeleri).

---

## 5. Görsel sistem (tokens) — `theme.py` yerine `tokens.py` + QSS şablonu

| Token | Koyu tema | Açık tema | Kullanım |
|---|---|---|---|
| `bg.base` | #0F1418 | #F5F7F8 | pencere |
| `bg.panel` | #161C21 | #FFFFFF | dock/panel |
| `bg.raised` | #1E262C | #EEF2F4 | kart, giriş alanı |
| `border` | #2E3A42 | #CBD5DB | kenarlık |
| `fg.primary` | #E8EEF1 | #14202A | metin (kontrast ≥ 12:1) |
| `fg.secondary` | #A9B8C1 | #4A5C68 | ikincil (≥ 6:1) |
| `fg.muted` | #7F909B | #6B7C87 | ipucu (≥ 4.5:1 — #78909a yerine) |
| `accent` | #2BC7CF | #0E8A93 | aktif adım, seçim |
| `ok` | #4DC56A | #1F8A3B | tamam |
| `warn` | #EFBC3F | #9A6B00 | uyarı |
| `error` | #EF6257 | #B3261E | hata/ret |
| `known` | #65D37C | — | overlay bilinen |
| `candidate` | #F3C348 | — | overlay aday |
| `rejected` | #EF6A60 | — | overlay ret |
| `neo` / `pha` | #FF9F43 / #FF4D4D | — | sınıf renkleri |
| `comet` / `kbo` / `planet` | #6EE7B7 / #C084FC / #FDE68A | — | sınıf renkleri |

Tipografi: temel 13 px (kullanıcı 12–16 arası ölçekleyebilir), başlık 15/18/22 px, tablo 12 px, caption 11 px **minimum**; font yığını `"Segoe UI", "Inter", "Noto Sans", system-ui`. Boşluk ölçeği 4/8/12/16/24 px. Köşe 6 px. Odak halkası 2 px `accent`.

QSS örneği (ölçekli değer üretimi):

```python
def build_qss(t: Tokens, base_px: int = 13) -> str:
    return f"""
    * {{ font-family: "Segoe UI", "Inter", "Noto Sans", system-ui; font-size: {base_px}px; color: {t.fg_primary}; }}
    QMainWindow, QWidget#appRoot {{ background: {t.bg_base}; }}
    QDockWidget::title {{ background: {t.bg_panel}; padding: 6px 10px; font-weight: 600; }}
    QPushButton {{ min-height: {int(base_px*2.3)}px; padding: 4px 12px; border-radius: 6px; background: {t.bg_raised}; border: 1px solid {t.border}; }}
    QPushButton:focus {{ border: 2px solid {t.accent}; }}
    QPushButton[role="primary"] {{ background: {t.accent}; color: {t.bg_base}; font-weight: 600; }}
    QLabel[tone="muted"] {{ color: {t.fg_muted}; }}
    QLabel[badge="ok"] {{ background: {t.ok}22; color: {t.ok}; border: 1px solid {t.ok}66; border-radius: 4px; padding: 3px 8px; }}
    ...
    """
```

---

## 6. Klavye ve verimlilik
`Ctrl+O` klasör · `Ctrl+R` yeniden çalıştır · `Ctrl+E` dışa aktar · `←/→` kare · `Home/End` ilk/son · `Space` blink · `B` blink aç/kapa · `D` fark modu · `F` sığdır · `1` 1:1 · `T` hedef seç modu · `C` karşılaştırma ekle · `A/R/N` kabul/ret/takip · `Ctrl+Z` inceleme kararını geri al · `Ctrl+L` günlük · `F1` bağlamsal yardım (aktif adımın kaide dokümanı bölümüne bağlantı).

---

## 7. Geçiş planı (mevcut koddan)
1. `theme.py` → `tokens.py` + `build_qss`; sabit piksel çağrılarını kaldır (`grep -n "setFixedWidth\|setMaximumWidth\|setFixedHeight"`).
2. `MainWindow`'u dock tabanlı iskelete çevir; `AsteroidWorkspace`/`ExoplanetWorkspace` içindeki panel oluşturucularını **adım widget'larına** böl (`steps/observation.py`, `steps/calibration.py`, …); ortak olanlar (`CalibrationStep`, `AstrometryStep`, `ImageViewer`, `LogDock`, `ProfileDialog`) paylaşılır.
3. Tüm metinleri `self.tr()` ile sar; `lupdate` ile `.ts` üret.
4. `QSettings` ile: son klasör, profil, dock yerleşimi, dil, ölçek, son seçilen yıldızlar (klasör hash'ine bağlı).
5. GUI testleri: `pytest-qt` ile adım geçişleri, boş durum, hata kartı; offscreen ekran görüntüsü regresyonu (`qa_preview/grab_ui.py` benzeri, her adım için).
