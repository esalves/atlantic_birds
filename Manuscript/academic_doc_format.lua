-- academic_doc_format.lua
-- Creates a dedicated Cover Page and Abstract Page for Word (.docx) documents,
-- with clean unnumbered title/abstract headers and OpenXML page breaks.
-- Preserves standard HTML rendering untouched.

function Pandoc(doc)
  if not quarto.doc.is_format("docx") then
    return doc
  end

  local meta = doc.meta
  local blocks = {}

  local function to_str(x)
    if not x then return "" end
    return pandoc.utils.stringify(x)
  end

  local pagebreak = pandoc.RawBlock('openxml', '<w:p><w:r><w:br w:type="page"/></w:r></w:p>')

  -- =========================================================================
  -- PAGE 1: COVER PAGE
  -- =========================================================================
  
  -- Title (unnumbered so that Pandoc doesn't prefix it with '1.')
  local title_text = to_str(meta.title)
  if title_text ~= "" then
    table.insert(blocks, pandoc.Header(1, { pandoc.Str(title_text) }, pandoc.Attr("", {"unnumbered"})))
  end

  -- Authors with Affiliation superscripts
  local author_inlines = {}
  local orcid_list = {}
  local corresponding_info = nil

  -- Quarto's manuscript engine organizes authors into meta['by-author']
  local authors = meta['by-author'] or meta.author or meta.authors

  if authors and type(authors) == "table" then
    for i, auth in ipairs(authors) do
      local name = ""
      if auth.name then
        if type(auth.name) == "table" and auth.name.literal then
          name = to_str(auth.name.literal)
        else
          name = to_str(auth.name)
        end
      end
      if name == "" then
        name = "Author " .. tostring(i)
      end

      -- Affiliation numbers
      local affil_num = tostring(i)
      local is_corr = false
      if auth.attributes and auth.attributes.corresponding then
        is_corr = true
      elseif auth.corresponding and (tostring(auth.corresponding) == "true" or auth.corresponding == true) then
        is_corr = true
      end

      local corr_marker = ""
      if is_corr then
        corr_marker = "*"
        local email = to_str(auth.email)
        if email ~= "" then
          corresponding_info = name .. " (" .. email .. ")"
        else
          corresponding_info = name
        end
      end

      local sup_str = affil_num .. corr_marker
      table.insert(author_inlines, pandoc.Str(name))
      table.insert(author_inlines, pandoc.Superscript({ pandoc.Str(sup_str) }))

      if i < #authors then
        table.insert(author_inlines, pandoc.Str(", "))
      end

      local orcid = to_str(auth.orcid)
      if orcid ~= "" then
        table.insert(orcid_list, name .. ": " .. orcid)
      end
    end
  end

  if #author_inlines > 0 then
    table.insert(blocks, pandoc.Para(author_inlines))
  end

  -- Spacing
  table.insert(blocks, pandoc.Para({ pandoc.Str("") }))

  -- Affiliations list
  local affils = meta['by-affiliation'] or meta.affiliations
  if affils and type(affils) == "table" then
    for i, aff in ipairs(affils) do
      local aff_name = ""
      if type(aff) == "table" and aff.name then
        aff_name = to_str(aff.name)
      else
        aff_name = to_str(aff)
      end
      if aff_name ~= "" then
        local aff_p = {
          pandoc.Superscript({ pandoc.Str(tostring(i)) }),
          pandoc.Space(),
          pandoc.Str(aff_name)
        }
        table.insert(blocks, pandoc.Para(aff_p))
      end
    end
  end

  -- Spacing
  table.insert(blocks, pandoc.Para({ pandoc.Str("") }))

  -- ORCIDs
  if #orcid_list > 0 then
    local orc_p = { pandoc.Strong({ pandoc.Str("ORCID: ") }) }
    for i, orc in ipairs(orcid_list) do
      table.insert(orc_p, pandoc.Str(orc))
      if i < #orcid_list then
        table.insert(orc_p, pandoc.Str(" | "))
      end
    end
    table.insert(blocks, pandoc.Para(orc_p))
  end

  -- Corresponding author
  if corresponding_info then
    local corr_p = {
      pandoc.Strong({ pandoc.Str("* Correspondence: ") }),
      pandoc.Str(corresponding_info)
    }
    table.insert(blocks, pandoc.Para(corr_p))
  end

  -- Short running title
  local running_p = {
    pandoc.Strong({ pandoc.Str("Running title: ") }),
    pandoc.Emph({ pandoc.Str("Morphological change in Atlantic Forest birds") })
  }
  table.insert(blocks, pandoc.Para(running_p))

  -- Page break at end of Cover Page
  table.insert(blocks, pagebreak)

  -- =========================================================================
  -- PAGE 2: ABSTRACT & KEYWORDS
  -- =========================================================================
  
  table.insert(blocks, pandoc.Header(1, { pandoc.Str("Abstract") }, pandoc.Attr("", {"unnumbered"})))

  if meta.abstract then
    if type(meta.abstract) == "table" and meta.abstract.t == "MetaBlocks" then
      for _, blk in ipairs(meta.abstract) do
        table.insert(blocks, blk)
      end
    else
      table.insert(blocks, pandoc.Para({ pandoc.Str(to_str(meta.abstract)) }))
    end
  end

  if meta.keywords then
    local kws = {}
    if type(meta.keywords) == "table" then
      for _, kw in ipairs(meta.keywords) do
        table.insert(kws, to_str(kw))
      end
    else
      table.insert(kws, to_str(meta.keywords))
    end
    local kw_p = {
      pandoc.Strong({ pandoc.Str("Keywords: ") }),
      pandoc.Str(table.concat(kws, ", "))
    }
    table.insert(blocks, pandoc.Para(kw_p))
  end

  if meta["plain-language-summary"] then
    table.insert(blocks, pandoc.Header(2, { pandoc.Str("Plain Language Summary") }, pandoc.Attr("", {"unnumbered"})))
    table.insert(blocks, pandoc.Para({ pandoc.Str(to_str(meta["plain-language-summary"])) }))
  end

  -- Page break at end of Abstract Page
  table.insert(blocks, pagebreak)

  -- =========================================================================
  -- PAGE 3+: MAIN TEXT
  -- =========================================================================
  for _, blk in ipairs(doc.blocks) do
    table.insert(blocks, blk)
  end

  -- Clear metadata so Pandoc does not emit duplicate unformatted headers
  meta.title = nil
  meta.author = nil
  meta.authors = nil
  meta['by-author'] = nil
  meta.date = nil
  meta.abstract = nil

  return pandoc.Pandoc(blocks, meta)
end
